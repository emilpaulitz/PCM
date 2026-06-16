import os

projDir = (os.path.abspath(os.getcwd())
           .removesuffix('/Code/gecko/panAraEcPcm/code')
           .removesuffix('/Code/gecko/panAraEcPcm')
           .removesuffix('/Code/gecko')
           .removesuffix('/Code') + '/')
geckoprojDir = projDir + 'Code/gecko/panAraEcPcm/'

# get accession from the files in projDir + 'Data/analysis/accSpecPCM/04_models/'
try:
    accs = [fname.removesuffix('.mat') 
            for fname in os.listdir(projDir + 'Data/analysis/accSpecPCM/04_models/')
            if fname.endswith('.mat')]

    for acc in accs:
        acc = acc.replace("-", "_")
        # put the adapter file into a folder named after the accession name, in the adapters folder
        # check if folder for adapter file is there; if not:
        if not os.path.exists(f'{geckoprojDir}adapters/{acc}'):
            os.makedirs(f'{geckoprojDir}adapters/{acc}')

        # copy the adapter_template and replace the string '$$acc$$' with the
        # accession name
        with open(f'{geckoprojDir}adapter_template.m', 'r') as template:
            template_str = template.read()
            new_str = template_str.replace('$$acc$$', acc).replace('$$projDir$$', projDir)
            with open(f'{geckoprojDir}adapters/{acc}/panAraEcPcmAdapter.m', 'w') as new_file:
                new_file.write(new_str)
except FileNotFoundError as e:
    print(e)
    print('This error can likely be fixed by running this script from the PCM folder.')

