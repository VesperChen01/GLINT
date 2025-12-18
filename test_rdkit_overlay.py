
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

def test_rdkit_overlay():
    mol = Chem.MolFromSmiles("c1ccccc1C(=O)NCC")
    Chem.AllChem.Compute2DCoords(mol)
    
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    opts = drawer.drawOptions()
    
    print("Step 1: Drawing Molecule...")
    drawer.DrawMolecule(mol)
    
    # 获取原子坐标
    atom_coords = []
    for i in range(mol.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords.append((pt.x, pt.y))
        print(f"Atom {i} coords: {pt.x}, {pt.y}")

    print("Step 2: Drawing Overlay Line...")
    drawer.SetColour((1.0, 0.0, 0.0, 1.0)) # Red
    drawer.SetLineWidth(5)
    # 画一条从左上到右下的对角线，跨越原子
    drawer.DrawLine(Point2D(10, 10), Point2D(590, 590))
    
    print("Step 3: Drawing Big Circle...")
    drawer.DrawEllipse(Point2D(250, 250), Point2D(350, 350))
    
    drawer.FinishDrawing()
    with open("test_overlay.png", "wb") as f:
        f.write(drawer.GetDrawingText())
    print("Saved test_overlay.png")

if __name__ == "__main__":
    test_rdkit_overlay()
