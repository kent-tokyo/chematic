use super::*;

#[test]
fn query_writer_is_idempotent_across_supported_grammar() {
    let cases = [
        "[H]",
        "[2H+]",
        "[13#6]",
        "[C,N,O]",
        "[!C]",
        "[$(C=O)]",
        "[N;H2,H1;!$(NC=O)]",
        "[C;D2;H1;+0]",
        "[C^3]",
        "[C;D{1-3}]",
        "[C;r{-6}]",
        "[C;R{2-}]",
        "[C;X2;v4;h1;x1;z1;Z1;k3]",
        "[C@H](F)(Cl)Br",
        "[C@@H](F)(Cl)Br",
        "C-C=C#N",
        "C$C",
        "c:c",
        "C~N",
        "C@C",
        "C!-N",
        "C-,=N",
        "C-&!@N",
        "F/C=C/F",
        "N->[Cu]",
        "[Cu]<-N",
        "C1CCCCC1",
        "C%12CCCCC%12",
        "C%(123)CCCCC%(123)",
        "C1(CC)CCC1",
        "C.C",
        "[si][as][se][te]",
        "[K,V,Y,W,U]",
    ];
    for input in cases {
        let once = smarts_to_smarts(input).unwrap_or_else(|e| panic!("{input}: {e}"));
        let twice = smarts_to_smarts(&once).unwrap_or_else(|e| panic!("{input} -> {once}: {e}"));
        assert_eq!(twice, once, "{input}");
    }
}

#[test]
fn malformed_and_unsupported_queries_are_rejected() {
    for input in [
        "C-", "-C", "C(", "C)", "C1CC", "[C", "[C: ]", "[C;D{1}]", "[C;D{-}]", "[C@TH1]", "C< C",
        "C->1CC1",
    ] {
        assert!(smarts_to_smarts(input).is_err(), "{input}");
    }
}

#[test]
fn reaction_writer_handles_agents_groups_and_dative_arrows() {
    let cases = [
        "C>>C",
        "C>O>N",
        "(C.O)>N>(C.N)",
        "N->[Cu]>O>N->[Cu]",
        "[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]",
    ];
    for input in cases {
        let once = reaction_to_smarts(input).unwrap_or_else(|e| panic!("{input}: {e}"));
        let twice = reaction_to_smarts(&once).unwrap_or_else(|e| panic!("{input} -> {once}: {e}"));
        assert_eq!(twice, once, "{input}");
    }

    for input in ["C>C", "C>>>C", "(C.O>C>C", "(C.O)N>C>C"] {
        assert!(reaction_to_smarts(input).is_err(), "{input}");
    }
}
