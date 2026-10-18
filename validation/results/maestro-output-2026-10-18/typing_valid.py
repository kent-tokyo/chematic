from chematic import Mol, MaeExport, write_mae_checked, v2000_to_mae_checked
def graph(mol: Mol) -> MaeExport:
    return write_mae_checked(mol, [[0., 0., 1.]], title="one", allow_metadata_loss=True)
def record(text: str) -> str:
    return v2000_to_mae_checked(text, allow_metadata_loss=True).text
