import json, sys
p = sys.argv[1]
nb = json.load(open(p, encoding='utf-8'))
n_err = 0
for i, c in enumerate(nb['cells']):
    if c['cell_type'] != 'code':
        continue
    for o in c.get('outputs', []):
        if o.get('output_type') == 'error':
            n_err += 1
            print(f"cell {i}: {o.get('ename')}: {o.get('evalue')}")
print("code cells:", sum(1 for c in nb['cells'] if c['cell_type']=='code'), "errors:", n_err)
