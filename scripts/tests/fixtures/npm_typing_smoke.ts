import init, { parse_smiles } from "@kent-tokyo/chematic";

export async function smoke(wasm: Uint8Array): Promise<string> {
  await init({ module_or_path: wasm });
  const molecule = parse_smiles("c1ccccc1");
  const formula: string = molecule.formula();
  molecule.free();
  return formula;
}
