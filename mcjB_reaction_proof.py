!pip install -q rdkit

import rdkit
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdChemReactions, Draw

# microcin J25 McjA precursor peptide MIKHFHFNKLSSGKKNNVPSPAKGVIQIKKSASQLTKGGAGHVPEYFVGIGTPISFYG
substrate = "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)CNC(=O)[C@H](C)NC(=O)CNC(=O)CNC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CO)NC(=O)[C@H](C)NC(=O)[C@H](CO)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@H](CCC(N)=O)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@@H]1CCCN1C(=O)[C@H](CO)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@@H](N)CCSC)[C@@H](C)CC)C(C)C)C(C)C)[C@@H](C)CC)[C@@H](C)CC)[C@@H](C)O)C(C)C)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O"

# reaction SMARTS for McjB, from MITE0000004
reaction_smarts = "[#6:1]-[#6:2](-[#6:3])-[#6:4]-[#6@H:5](-[#7:6]-[#6:7](=[O:8])-[#6@H:9](-[#6:10]-[#6:11]-[#6:12](-[#7:13])=[O:14])-[#7:15]-[#6:16](=[O:17])-[#6@H:18](-[#6:19]-[#8:20])-[#7:21]-[#6:22](=[O:23])-[#6@H:24](-[#6:25])-[#7:26]-[#6:27](=[O:28])-[#6@H:29](-[#6:30]-[#8:31])-[#7:32]-[#6:33](=[O:34])-[#6@@H:35](-[#7:36])-[#6:37]-[#6:38]-[#6:39]-[#6:40]-[#7:41])-[#6:42](=[O:43])-[#7:44]-[#6@@H:45](-[#6@@H:62](-[#6:63])-[#8:64])-[#6:46](=[O:47])-[#7:48]-[#6@@H:49](-[#6:50]-[#6:51]-[#6:52]-[#6:53]-[#7:54])-[#6:55](=[O:56])-[#7:57]-[#6:58]-[#6:59](-[#7:61])=[O:60]>>[#7:57]-[#6:58]-[#6:59](-[#7:61])=[O:60]"

# intermediate core peptide GGAGHVPEYFVGIGTPISFYG
expected_product = "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)CNC(=O)[C@H](C)NC(=O)CNC(=O)CN)C(C)C)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O"

substrate = Chem.MolFromSmiles(substrate)
substrate = Chem.RemoveHs(substrate)
expected_product = Chem.MolFromSmiles(expected_product)
reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)

products = reaction.RunReactants([substrate])
for product in products:
    Chem.SanitizeMol(product[0])

if (
    len(set([Chem.MolToSmiles(product[0]) for product in products])) == 1 and
    Chem.MolToSmiles(products[0][0]) == Chem.MolToSmiles(expected_product)
):
    print("The reaction product matches the expected product.")
    print(f"{Chem.MolToSmiles(products[0][0])}")
