!pip install -q rdkit

import rdkit
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdChemReactions, Draw

# intermediate core peptide GGAGHVPEYFVGIGTPISFYG
substrate = "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)CNC(=O)[C@H](C)NC(=O)CNC(=O)CN)C(C)C)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O"

# reaction SMARTS for McjC, from MITE0000001
reaction_smarts = "[#6:55]-[#6:54](-[#7:56]-[#6:57](=[O:58])-[#6:59]-[#7:60]-[#6:61](=[O:62])-[#6:63]-[#7:64])-[#6:52](=[O:53])-[#7:51]-[#6:50]-[#6:48](=[O:49])-[#7:47]-[#6@@H:40](-[#6:41]-[c:42]1[c:43][n:44][c:45][n:46]1)-[#6:38](=[O:39])-[#7:37]-[#6@@H:36](-[#6:65])-[#6:34](=[O:35])-[#7:33]-[#6:32]-[#6:30](=[O:31])-[#7:29]-[#6@@H:23](-[#6:24]-[#6:25]-[#6:26](-[#8:28])=[O:27])-[#6:21](=[O:22])-[#7:20]-[#6@@H:18](-[#6:19])-[#6:16](=[O:17])-[#7:15]-[#6@@H:13](-[#6:14])-[#6:11](=[O:12])-[#7:10]-[#6:9]-[#6:7](=[O:8])-[#7:6]-[#6:5]-[#6:3](=[O:4])-[#7:2]-[#6:1]-[#6:66](=[O:67])-[#7:68]-[#6:69]-[#6:70](=[O:71])-[#7:72]-[#6:73]-[#6:74](=[O:75])-[#7:76]-[#6@@H:77](-[#6:111])-[#6:78](=[O:79])-[#7:80]-[#6@@H:81](-[#6:110])-[#6:82](=[O:83])-[#7:84]-[#6@@H:85](-[#6:86])-[#6:87](=[O:88])-[#7:89]-[#6@@H:90](-[#6:91])-[#6:92](=[O:93])-[#7:94]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:103](=[O:104])-[#7:105]-[#6:106]-[#6:107](-[#8:109])=[O:108]>>[#6:111]-[#6@H:77](-[#7:76]-[#6:74](=[O:75])-[#6:73]-[#7:72]-[#6:70](=[O:71])-[#6:69]-[#7:68]-[#6:66](=[O:67])-[#6:1]-[#7:2]-[#6:3](=[O:4])-[#6:5]-[#7:6]-[#6:7](=[O:8])-[#6:9]-[#7:10]-[#6:11](=[O:12])-[#6@H:13](-[#6:14])-[#7:15]-[#6:16](=[O:17])-[#6@H:18](-[#6:19])-[#7:20]-[#6:21](=[O:22])-[#6@@H:23]-1-[#6:24]-[#6:25]-[#6:26](=[O:27])-[#7:64]-[#6:63]-[#6:61](=[O:62])-[#7:60]-[#6:59]-[#6:57](=[O:58])-[#7:56]-[#6:54](-[#6:55])-[#6:52](=[O:53])-[#7:51]-[#6:50]-[#6:48](=[O:49])-[#7:47]-[#6@@H:40](-[#6:41]-[c:42]2[c:43][n:44][c:45][n:46]2)-[#6:38](=[O:39])-[#7:37]-[#6@@H:36](-[#6:65])-[#6:34](=[O:35])-[#7:33]-[#6:32]-[#6:30](=[O:31])-[#7:29]-1)-[#6:78](=[O:79])-[#7:80]-[#6@@H:81](-[#6:110])-[#6:82](=[O:83])-[#7:84]-[#6@@H:85](-[#6:86])-[#6:87](=[O:88])-[#7:89]-[#6@@H:90](-[#6:91])-[#6:92](=[O:93])-[#7:94]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:103](=[O:104])-[#7:105]-[#6:106]-[#6:107](-[#8:109])=[O:108]"

# mature microcin J25
expected_product = "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H]1CCC(=O)NCC(=O)NCC(=O)N[C@@H](C)C(=O)NCC(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](C(C)C)C(=O)N2CCC[C@H]2C(=O)N1)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O"

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
