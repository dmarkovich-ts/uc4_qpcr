# import asyncio

from csv import QUOTE_NONNUMERIC
from json import dump, load
from numpy import concatenate, unique
from math import log10
from pandas import DataFrame, merge, notna
from pprint import pprint
from time import time

from LuPySNB.API.MaterialLibrary import get_material, get_material_field_values
from LuPySNB.API.Plate import get_plate_container, \
    get_plates_from_plate_container
from LuPySNB.REST import RestClientHTTPX, RestConfig


def get_client():
    # Get client for connecting to Signals

    with open('.signals/auth.json') as fp:
        auth = load(fp)

    headers = {
        "x-api-key": auth['api_key'],
        "Content-Type": "application/vnd.api+json",
        "Accept": "application/vnd.api+json"
    }

    config = RestConfig(
        base_url=auth['base_url'],
        headers=headers, auto_expand=True, next_limit=15
    )
    return RestClientHTTPX(config)


def get_experiment(plate_container):
    # Get metadata of experiment that containes the plate container

    parents = [item for item in plate_container.included
               if item.attributes.type == 'experiment']
    if len(parents) == 1:
        return {'eid': parents[0].id, 'name': parents[0].attributes.name}
    else:
        return {'eid': '', 'name': ''}


def get_plate_container_name(plate_container, plate_container_eid):
    for item in plate_container.included:
        if item.id == plate_container_eid and item.type == 'entity':
            # TODO check that name is available in attributes
            return item.attributes.name
    return ''


def row2letter(row: int) -> str:
    """Map a positive integer to a base-26 letter sequence (A, B, ...,
    Z, AA, AB, ...).

    1 -> A
    26 -> Z
    27 -> AA

    """
    if row < 1:
        raise ValueError("n must be a positive integer (starting from 1)")

    letters = []
    while row > 0:
        row, remainder = divmod(row - 1, 26)
        letters.append(chr(ord('A') + remainder))
    return ''.join(reversed(letters))


def get_plates(client, plate_container_eid):
    # get plates from plate_container_eid as list of DataFrame

    plate_container = get_plate_container(client, plate_container_eid)
    plate_container_name = get_plate_container_name(plate_container,
                                                    plate_container_eid)

    experiment = get_experiment(plate_container)
    pprint(experiment)

    plate_list_response = \
        get_plates_from_plate_container(client, plate_container_eid)

    plates = []
    for plate in plate_list_response.data:
        print(f'>> Processing plate {plate.attributes.name}, '
              f'eid = {plate.attributes.eid}')

        tmp = {
            'name': plate.attributes.name,
            'attributes': {'n_rows': plate.attributes.numberOfRows,
                           'n_columns': plate.attributes.numberOfColumns},
            'data': None,
            'experiment': experiment,
            'plate_container': {
                'eid': plate_container_eid,
                'name': plate_container_name
            }
        }

        data = []  # list of dict
        for annotation_layer in plate.attributes.annotationLayers:
            # skip metadata for LU Plate Editor Signals app
            if not annotation_layer.name == 'plate_editor_data':
                print(f' > Processing annotation layer {annotation_layer.name}'
                      f' of type {annotation_layer.dataType}')
                if annotation_layer.dataType == 'ATTRIBUTE':
                    # attributes are global for entire plate, avoid duplicating
                    tmp["attributes"].update(
                        {annotation_layer.name: annotation_layer.value}
                    )
                elif annotation_layer.dataType == "LINK":
                    for well in annotation_layer.wells:
                        # if well and well.content:
                        well_content = True if well.content else False
                        data.append({
                            'wellId': well.wellId if well.wellId else None,
                            'column': annotation_layer.name,
                            'user': well.content.user if well_content else None,
                            'value': well.content.value if well_content else None
                        })
                    # End for well
                else:
                    print(f'TODO unknown type {annotation_layer.dataType}')
        # End for annotation_layer

        # unfilled wells will be omitted in annotation layers, fix by
        # creating a DF with proper zero-filled well identifiers and merge
        layout = []
        zfill_len = 1 + int(log10(tmp['attributes']['n_columns']))
        for row in range(tmp['attributes']['n_rows']):
            for col in range(tmp['attributes']['n_columns']):
                layout.append({
                    'wellId': f'{row2letter(row + 1) + str(col + 1)}',
                    'well_id': f'{row2letter(row + 1) + str(col + 1).zfill(zfill_len)}'
                })
        layout = DataFrame.from_records(layout)  # merge later woth data

        tmp["data"] = merge(left=DataFrame.from_records(data), right=layout,
                            how='outer', on='wellId')
        tmp["data"].sort_values(by='well_id', inplace=True)

        plates.append(tmp)
        # End for plate
    return plates


rules = {
    "resolve": {
        'Biosample Containers': {
            'tokenized_expression': [
                'Biosample Container Name', ' : ', 'Animal Batch.Species'
            ],
            'derived_tokenized_expression': ['Matrix']
        },
        'Cell Sample Containers': {
            'tokenized_expression': [
                'Cell Sample Container Name'
            ],
            'derived_tokenized_expression': ['Barcode']
        },
        'Lu Primers': {
            'tokenized_expression': ['Fluorophore (with quenchers)'],
            'derived_tokenized_expression': ['Gene Name']
        }
    },
    "derived_columns_mapping": {
        "Sample": 'Biological Set Name',
        'Channel 1': 'Channel 1 Gene Name',
        'Channel 2': 'Channel 2 Gene Name',
        'Channel 3': 'Channel 3 Gene Name',
        'Channel 4': 'Channel 4 Gene Name',
        'Channel 5': 'Channel 5 Gene Name'
    }
}


def resolve_tokenized_expression(tokenized_expression, material_asset_response,
                                 client):
    # Resolves asset link to value as described by tokenized expression

    keys = material_asset_response.data.attributes.fields.keys()

    resolved_value = []
    for tkn in tokenized_expression:
        if tkn in keys or tkn.split('.')[0] in keys:
            resolved_value.append(get_material_field_values(
                client, material_asset_response, tkn
            ))
        else:
            # (parent) tkn NOT found in asset fields are delimiters
            # that should be just pasted into the resulting string
            resolved_value.append(tkn)

    return ''.join(resolved_value)


def resolve_unique_assets(plates, rules, client):
    # Gets unique assets from `plates`, resolves them according to
    # `rules`, returns a DataFrame of resolved values and derived
    # columns / resolved values

    unique_assets = unique(concatenate([
        plate["data"]["value"].dropna().unique() for plate in plates
    ])).tolist()
    print(f'>> Found n = {len(unique_assets)} unique assets to resolve')

    assets_df = []
    # TODO consider using `_async` LuPySNB methods or parallelization
    for asset in unique_assets:
        print(f'>> Resolving asset = {asset}')
        # start = time()
        material_asset_response = get_material(client, asset)
        lib = material_asset_response.data.attributes.library

        if lib not in rules.keys():
            print(f'>> Unknown lib = {lib} for asset = {asset}')

        resolved_value = resolve_tokenized_expression(
            rules[lib]['tokenized_expression'],
            material_asset_response, client
        )
        # derived_column_postfix = rules[lib]['derived']['column_postfix']
        derived_column_resolved_value = resolve_tokenized_expression(
            rules[lib]['derived_tokenized_expression'],
            material_asset_response, client
        )
        # print(f'   Done in {time() - start} seconds')  # [0.3, 1] seconds
        assets_df.append({
            'value': asset,
            'resolved_value': resolved_value,
            # 'derived_column_postfix': derived_column_postfix,
            'derived_column_resolved_value': derived_column_resolved_value
        })
    # End for asset

    return DataFrame.from_records(assets_df)


def resolve_plates(plates, resolved_unique_assets, derived_columns_mapping):
    # Modifies field 'data' of each plate, resolving values

    mpng = DataFrame()  # maps derived resolved values to columns
    mpng['column'] = derived_columns_mapping.keys()
    mpng['derived_column'] = derived_columns_mapping.values()

    for plate in plates:
        tmp = merge(left=plate['data'], right=resolved_unique_assets,
                    how='left', on='value')
        tmp = merge(left=tmp, right=mpng, how='left', on='column')
        # tmp['derived_column'] = tmp['column'] + tmp['derived_column_postfix']
        # tmp.head()
        # tmp.info()

        # Keep sample asset until Replicate column is constructed: if
        # two different links result in identical resolved value,
        # samples should still be assigned different Replicate indexes
        sample_asset = tmp.pivot(columns='column', index='well_id',
                                 values='value')[['Sample']]
        sample_asset.rename(columns={'Sample': 'sample_asset'}, inplace=True)
        # sample_asset.head()
        # sample_asset.info()

        resolved = tmp.pivot(columns='column', index='well_id',
                             values='resolved_value')
        resolved.dropna(axis='columns', how='all', inplace=True)
        # resolved.head()
        # resolved.info()

        derived = tmp.pivot(columns='derived_column', index='well_id',
                            values='derived_column_resolved_value')
        derived.dropna(axis='columns', how='all', inplace=True)
        # derived.head()
        # derived.info()

        full = merge(left=resolved, right=derived, how='left',
                     left_index=True, right_index=True)
        full = merge(left=full, right=sample_asset, how='left',
                     left_index=True, right_index=True)
        full.reset_index(inplace=True)
        # full.head()
        # full.info()

        # Wells belong to same Replicate index if they have identical
        # values in columns `repl_clmns`:
        repl_clmns = full.columns.tolist()
        for clmn in ['Sample', 'well_id']:
            repl_clmns.remove(clmn)
        # Get all unique combinations of `repl_clmns` and enumerate them from 1
        repl_df = full[notna(full.sample_asset)]\
            .groupby(repl_clmns, dropna=False)\
            .size().reset_index().rename(columns={0: 'Replicate'})
        repl_df['Replicate'] = 1 + repl_df.index
        # repl_df.head()
        # repl_df.info()

        full = merge(left=full, right=repl_df, how='left', on=repl_clmns)
        full.Replicate = full.Replicate.astype('Int64')
        # full.head()
        # full.info()

        plate['data'] = full
    # End for plate
    return None


def get_asset_id(asset):
    # Returns id part of asset string `<asset>:<id>`
    return asset.split(":")[1]


def get_plrn_df(df, mpng):

    # Creates DataFrame corresponding to what section `Plate Data` of PLRN files expects

    # plate_id = df.iloc[0]['Plate']
    # instrument_protocol = df.iloc[0]['Instrument Protocol']

    clmns = df.columns.tolist()

    plrn_df = DataFrame()
    plrn_df['Well'] = df['well_id']

    for i in range(1, 6, 1):
        clmn = f'Channel {i}'
        plrn_df[f'Ch{i} Dye'] = df[clmn] if clmn in clmns else None

    plrn_df['FRET'] = None

    plrn_df['Sample Type'] = 'Unknown'
    plrn_df['Sample Name'] = df['Sample']

    for i in range(1, 6, 1):
        clmn = mpng[f'Channel {i}']
        plrn_df[f'Ch{i} Target Name'] = df[clmn] if clmn in clmns else None
    plrn_df['FRET Target Name'] = None

    plrn_df['Biological Set Name'] = df[mpng['Sample']]
    plrn_df['Replicate'] = df['Replicate']

    for i in range(1, 6, 1):
        plrn_df[f'Ch{i} Quantity'] = None
    plrn_df['FRET Quantity'] = None

    plrn_df['Well Note'] = None
    for i in range(1, 6, 1):
        plrn_df[f'Ch{i} Well Color'] = None
    plrn_df['FRET Well Color'] = None

    return plrn_df


def plrn(plate, plrn_df):
    experiment_eid = plate['experiment']['eid']
    container_name = plate['plate_container']['name']
    plate_name = plate['name']
    instrument_protocol = plate['attributes']['Instrument Protocol'] \
        if 'Instrument Protocol' in plate['attributes'] else ''
    return (
        f'Plate Header,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Filed,Data,,Instruction,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Version,1,,Do not modify this field.,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Plate Size,96,,Do not modify this field.,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Plate Type,BR Clear,,Allowed values (BR White,BR Clear),,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Scan Mode,All Channels,,"Allowed values (""SYBR/FAM Only"""," ""All Channels"""," ""FRET"")",,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Units,copy number,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Run ID,{experiment_eid},{container_name},{plate_name},,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Run Notes,{instrument_protocol},,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Run Protocol,RespProtocol.prcl,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Data File,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'TBD,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'Plate Data,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n'
        f'{plrn_df.to_csv(index=False, quoting=QUOTE_NONNUMERIC)}'
    )


if __name__ == '__main__':

    # Test examples from SNB:EXP-134
    # https://lundbeck-sandbox.signalsresearch.revvitycloud.eu/elements/entity/experiment:42c9df43-3da6-466d-82c9-9dc3579de1d1?focus=experiment%3Ab6f190c9-ea2a-4ec2-a0ee-54991ba8609c

    # One_Plate_example
    plate_container_eid = 'plateContainer:ed19a9aa-6ae6-4c93-b170-45bf6755bb75'
    # Two_Plate_example
    # plate_container_eid = 'plateContainer:e0c0e8c0-d057-49e7-8953-1cc44dca286e'
    # Cell samples plate
    # plate_container_eid = 'plateContainer:fd77e48f-bbb2-4612-a1fb-96f9b50ff52f'
    # Swiss Cheese plate
    # plate_container_eid = 'plateContainer:0715e4d5-aa28-448e-96e2-28e1d45c0b54'

    client = get_client()

    # 1. read plates and save unresolved results

    plates = get_plates(client, plate_container_eid)
    pprint(plates)

    for plate in plates:
    # TODO labels on saved files for experiment_name and instrument_type
        fname = f'data/forward/experiment__{plate["experiment"]["name"]}__'\
            f'plate_container__{plate["plate_container"]["name"]}__'\
            f'plate__{plate["name"]}__unresolved.json'
        print(f'>> Saving unresolved plate {fname}')
        cntnt = plate.copy()
        cntnt['data'] = cntnt['data'].to_dict(orient='list')
        with open(fname, 'w') as fp:
            dump(cntnt, fp, indent=2)

    # 2. check if individual plates are for Bio-Rad, if yes convert to PLRN

    for plate in plates:
        if 'Instrument Type' in plate['attributes'] and \
           plate['attributes']['Instrument Type'] == 'Bio Rad / CFX96 / 5 Channel':
            print(f'>> Performing PLRN conversion for plate {plate["name"]} with rules:')
            pprint(rules)

            # TODO resolving for multiple plates at the same time is
            # much more efficient, as they would likely share many
            # unique assets. The pipelines will not work on a single
            # plate at a time, for implementation compatibility plate
            # objects are passed as single-element lists
            resolved_unique_assets = resolve_unique_assets([plate],
                                                           rules['resolve'], client)
            print(resolved_unique_assets)

            resolve_plates([plate], resolved_unique_assets,
                           derived_columns_mapping=rules['derived_columns_mapping'])
            print(' > Resolved plate:')
            pprint(plate)

            fname_base = f'data/forward/experiment__{plate["experiment"]["name"]}__'\
                f'plate_container__{plate["plate_container"]["name"]}__'\
                f'plate__{plate["name"]}__resolved'
            plate['data'].to_csv(path_or_buf=f'{fname_base}.csv',
                                 index=False, quoting=QUOTE_NONNUMERIC)

            plrn_df = get_plrn_df(plate['data'],
                                  rules['derived_columns_mapping'])
            print(' > Resolved PLRN:')
            print(plrn_df)

            with open(f'{fname_base}.plrn', 'w') as fl:
                fl.write(plrn(plate, plrn_df))
