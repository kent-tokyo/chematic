package io.chematic.validation;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.BufferedReader;
import java.io.StringReader;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import org.openscience.cdk.CDK;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public final class CdkAdapter {
  private static final ObjectMapper JSON = new ObjectMapper();
  private static final String[] SMARTS = {
      "[#6]", "[#7]", "[#8]", "c1ccccc1", "[C,c](=O)[O,N]", "[N;H1,H2,H3;+0]"
  };
  private static final SmilesParser PARSER =
      new SmilesParser(DefaultChemObjectBuilder.getInstance());
  private static final SmilesGenerator SMILES =
      new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.Stereo | SmiFlavor.UseAromaticSymbols);
  private static final Map<String, SmartsPattern> QUERIES = new LinkedHashMap<>();
  private static final Map<String, Function<IAtomContainer, Object>> OPERATIONS =
      new LinkedHashMap<>();

  static {
    for (String query : SMARTS) QUERIES.put(query, SmartsPattern.create(query));
    OPERATIONS.put("canonical_stable", mol -> canonicalStable(mol));
    OPERATIONS.put("formula", CdkAdapter::formula);
    OPERATIONS.put("molecular_weight_centi_da", mol -> Math.round(descriptorDouble(new WeightDescriptor(), mol) * 100));
    OPERATIONS.put("tpsa_milli", mol -> Math.round(descriptorDouble(new TPSADescriptor(), mol) * 1000));
    OPERATIONS.put("hba", mol -> descriptorInteger(new HBondAcceptorCountDescriptor(), mol));
    OPERATIONS.put("hbd", mol -> descriptorInteger(new HBondDonorCountDescriptor(), mol));
    OPERATIONS.put("rings", mol -> Cycles.sssr(mol).numberOfCycles());
    OPERATIONS.put("mol_roundtrip", CdkAdapter::molRoundtrip);
    for (Map.Entry<String, SmartsPattern> entry : QUERIES.entrySet()) {
      OPERATIONS.put("smarts:" + entry.getKey(), mol -> entry.getValue().matches(mol));
    }
  }

  private CdkAdapter() {}

  private static IAtomContainer parse(String smiles) {
    try {
      return PARSER.parseSmiles(smiles);
    } catch (Exception exc) {
      throw new IllegalArgumentException(exc);
    }
  }

  private static String canonical(IAtomContainer mol) {
    try {
      return SMILES.create(mol);
    } catch (CDKException exc) {
      throw new IllegalStateException(exc);
    }
  }

  private static boolean canonicalStable(IAtomContainer mol) {
    String text = canonical(mol);
    return canonical(parse(text)).equals(text);
  }

  private static String formula(IAtomContainer mol) {
    String value = MolecularFormulaManipulator.getString(
        MolecularFormulaManipulator.getMolecularFormula(mol));
    if (value.startsWith("[") && value.contains("]")) {
      int close = value.indexOf(']');
      value = value.substring(1, close) + value.substring(close + 1);
    }
    return value;
  }

  private static double descriptorDouble(Object descriptor, IAtomContainer mol) {
    DescriptorValue value;
    if (descriptor instanceof WeightDescriptor d) value = d.calculate(mol);
    else if (descriptor instanceof TPSADescriptor d) value = d.calculate(mol);
    else throw new IllegalArgumentException("unsupported double descriptor");
    return ((DoubleResult) value.getValue()).doubleValue();
  }

  private static int descriptorInteger(Object descriptor, IAtomContainer mol) {
    DescriptorValue value;
    if (descriptor instanceof HBondAcceptorCountDescriptor d) value = d.calculate(mol);
    else if (descriptor instanceof HBondDonorCountDescriptor d) value = d.calculate(mol);
    else throw new IllegalArgumentException("unsupported integer descriptor");
    return ((IntegerResult) value.getValue()).intValue();
  }

  private static boolean molRoundtrip(IAtomContainer mol) {
    try {
      StringWriter buffer = new StringWriter();
      try (MDLV2000Writer writer = new MDLV2000Writer(buffer)) {
        writer.write(mol);
      }
      IAtomContainer reread = DefaultChemObjectBuilder.getInstance().newAtomContainer();
      try (MDLV2000Reader reader = new MDLV2000Reader(new StringReader(buffer.toString()))) {
        reread = reader.read(reread);
      }
      AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(reread);
      Aromaticity.cdkLegacy().apply(reread);
      return canonical(reread).equals(canonical(mol));
    } catch (Exception exc) {
      throw new IllegalStateException(exc);
    }
  }

  private static String sha256(byte[] value) throws Exception {
    return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(value));
  }

  private static Map<String, Object> operationResult(IAtomContainer mol, String operation) {
    Map<String, Object> result = new LinkedHashMap<>();
    try {
      result.put("status", "ok");
      result.put("value", OPERATIONS.get(operation).apply(mol));
    } catch (Exception exc) {
      result.put("status", "error");
      result.put("error", exc.getClass().getSimpleName() + ": " + exc.getMessage());
    }
    return result;
  }

  private static void emitResults(Path corpus) throws Exception {
    byte[] bytes = Files.readAllBytes(corpus);
    String corpusHash = sha256(bytes);
    try (BufferedReader reader = Files.newBufferedReader(corpus, StandardCharsets.UTF_8)) {
      String line;
      while ((line = reader.readLine()) != null) {
        if (line.isBlank()) continue;
        Map<String, Object> row = JSON.readValue(line, new TypeReference<>() {});
        Map<String, Object> record = new LinkedHashMap<>();
        record.put("schema_version", 1);
        record.put("engine", "cdk");
        record.put("engine_version", CDK.getVersion());
        record.put("source_commit", null);
        record.put("corpus_sha256", corpusHash);
        record.put("id", row.get("id"));
        record.put("smiles", row.get("smiles"));
        record.put("status", "ok");
        Map<String, Object> operations = new LinkedHashMap<>();
        try {
          String smiles = (String) row.get("smiles");
          parse(smiles);
          for (String operation : OPERATIONS.keySet()) {
            // CDK descriptors and writers are allowed to mutate atom-container
            // state. Accuracy rows must not depend on which operation ran first.
            operations.put(operation, operationResult(parse(smiles), operation));
          }
        } catch (Exception exc) {
          record.put("status", "parse_error");
          operations.put("parse", Map.of(
              "status", "parse_error", "error", exc.getClass().getSimpleName() + ": " + exc.getMessage()));
        }
        record.put("operations", operations);
        System.out.println(JSON.writeValueAsString(record));
      }
    }
  }

  private static Map<String, Object> benchmark(Map<String, Object> request) throws Exception {
    @SuppressWarnings("unchecked")
    List<String> smiles = (List<String>) request.get("smiles");
    String operation = (String) request.get("operation");
    String lane = (String) request.get("lane");
    int iterations = ((Number) request.getOrDefault("iterations", 1)).intValue();
    if (!operation.equals("parse") && !OPERATIONS.containsKey(operation)) {
      throw new IllegalArgumentException("unsupported operation: " + operation);
    }
    List<IAtomContainer> prepared = null;
    if (lane.equals("prepared")) {
      prepared = new ArrayList<>();
      for (String text : smiles) prepared.add(parse(text));
    } else if (!lane.equals("pipeline")) {
      throw new IllegalArgumentException("unsupported lane: " + lane);
    }
    List<Object> values = new ArrayList<>();
    int errors = 0;
    long start = System.nanoTime();
    for (int iteration = 0; iteration < iterations; iteration++) {
      for (int index = 0; index < smiles.size(); index++) {
        try {
          IAtomContainer mol = prepared == null ? parse(smiles.get(index)) : prepared.get(index);
          values.add(operation.equals("parse") ? true : OPERATIONS.get(operation).apply(mol));
        } catch (Exception exc) {
          errors++;
          values.add(List.of("error", exc.getClass().getSimpleName()));
        }
      }
    }
    long elapsed = System.nanoTime() - start;
    Map<String, Object> response = new LinkedHashMap<>();
    response.put("engine", "cdk");
    response.put("engine_version", CDK.getVersion());
    response.put("operation", operation);
    response.put("lane", lane);
    response.put("rows", smiles.size() * iterations);
    response.put("elapsed_ns", elapsed);
    response.put("errors", errors);
    response.put("digest", sha256(JSON.writeValueAsBytes(values)));
    return response;
  }

  private static void serve() throws Exception {
    try (BufferedReader reader = new BufferedReader(new java.io.InputStreamReader(System.in))) {
      String line;
      while ((line = reader.readLine()) != null) {
        Map<String, Object> response;
        try {
          Map<String, Object> request = JSON.readValue(line, new TypeReference<>() {});
          String command = (String) request.get("command");
          if ("metadata".equals(command)) {
            List<String> operations = new ArrayList<>();
            operations.add("parse");
            operations.addAll(OPERATIONS.keySet());
            response = Map.of(
                "engine", "cdk", "engine_version", CDK.getVersion(), "operations", operations);
          } else if ("benchmark".equals(command)) {
            response = benchmark(request);
          } else {
            throw new IllegalArgumentException("unknown command: " + command);
          }
        } catch (Exception exc) {
          response = Map.of("error", exc.getClass().getSimpleName() + ": " + exc.getMessage());
        }
        System.out.println(JSON.writeValueAsString(response));
        System.out.flush();
      }
    }
  }

  public static void main(String[] args) throws Exception {
    if (args.length == 1 && args[0].equals("--server")) {
      serve();
      return;
    }
    if (args.length == 2 && args[0].equals("--corpus")) {
      emitResults(Path.of(args[1]));
      return;
    }
    throw new IllegalArgumentException("usage: --server | --corpus PATH");
  }
}
