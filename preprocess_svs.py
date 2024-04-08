import argparse
import os
import json
from tqdm import tqdm

import numpy as np
import pandas as pd

OPENSLIDE_BIN_PATH = r'C:\Program Files\openslide\bin'
if hasattr(os, 'add_dll_directory'):
    # Windows
    with os.add_dll_directory(OPENSLIDE_BIN_PATH):
        import openslide
else:
    import openslide

def get_diag_and_treatm_per_patient(cd_json):
    cd_df = None
    with open(cd_json, 'r') as f:
        clinical_data = json.load(f)
        max_len = 0
        cd_per_pat = {}
        for i in range(len(clinical_data)):
            try:
                id_ = clinical_data[i]['diagnoses'][0]['submitter_id'].strip('_diagnosis')
                cd_per_pat[id_] = [clinical_data[i]['diagnoses'][0]['ajcc_pathologic_stage']]
                cd_per_pat[id_].append(clinical_data[i]['diagnoses'][0]['site_of_resection_or_biopsy'])
                #cd_per_pat[id_].append(clinical_data[i]['diagnoses'][0]['state'])
                #cd_per_pat[id_].append(len(clinical_data[i]['diagnoses']))
                cd_per_pat[id_].append(clinical_data[i]['project']['project_id'])
                cd_per_pat[id_].append(clinical_data[i]['case_id'])
                #if 'treatments' in clinical_data[i]['diagnoses'][0]:
                #    cd_per_pat[id_].extend([t['treatment_type'] for t in clinical_data[i]['diagnoses'][0]['treatments'] if 'treatment_type' in t])
                if len(cd_per_pat[id_]) > max_len:
                    max_len = len(cd_per_pat[id_])
            except KeyError as e:  # stage not reported
                if any('ajcc_pathologic_stage' in d for d in clinical_data[i]['diagnoses']):
                    print(id_, e)
                continue

        for id_ in cd_per_pat.keys():
            cd_per_pat[id_].extend(['' for j in range(max_len - len(cd_per_pat[id_]))])
        cd_df = pd.DataFrame(cd_per_pat, index=['stage', 'biopsy site', 'project_id', 'case_id']).T
    return cd_df


def extract_slide_info(p, cd_df, results_dir):
    '''
    requires `cd_df` to be set, from `get_diag_and_treatm_per_patient`
    extracts and stores slide image to `results_dir`
    returns meta data from svs path (patient-specific) and carcinoma stage from primary diagnosis
    '''
    id_ = '-'.join(p.split('\\')[-1].split('-')[:3])
    try:
        stage = cd_df.loc[id_, 'stage']
    except KeyError:
        return None

    slide = openslide.OpenSlide(p)
    image = slide.associated_images['thumbnail']
    image_dir_path = results_dir + '\\' +  cd_df.loc[id_, 'stage']
    os.makedirs(image_dir_path, exist_ok=True)
    image_path = image_dir_path + '\\' + id_ + '.png'
    image.save(image_path, 'PNG')

    size = np.array(image).shape

    slide_info_ser = pd.Series({
        'id': id_,
        'size': size,
        'path': p,
        'image_path': image_path,
        **dict(slide.properties),
        **cd_df.loc[id_].to_dict()
    })
    return slide_info_ser

def extract_image_and_stage_from_svs(root_dir, cd_json, results_dir):
    '''
    Extracting images and carcinoma stage from .svs files and cohort clinical data json
    '''
    # clinical data of cohort
    cd_df = get_diag_and_treatm_per_patient(cd_json)

    # .svs paths
    paths = []
    for path, dirs, files in os.walk(root_dir):
        svs_files = [path + '\\' + f for f in files if f.endswith('.svs')]
        paths.extend(svs_files)

    # extracts clinical data and WSI meta-data to `slides_info`
    # dumps WSI from .svs to `results_dir`
    slides_info = []
    for p in tqdm(paths):
        slide_info = extract_slide_info(p, cd_df, results_dir)
        slides_info.append(slide_info)

    slides_info = pd.concat(slides_info, axis=1).T
    slides_info.to_csv(results_dir + '\\slides_info.csv')
    return slides_info


def main(args):
    print('Extraction started ...')
    slides_info = extract_image_and_stage_from_svs(args.root, args.clinical_data_json, args.results)
    print('Extraction done!')
    print(slides_info)


if __name__ == "__main__":
    def_results_dir = r'C:\Users\DELL\Documents\Github\carcinoma_stage_detection_wsi\data'
    def_root_dir = r'C:\Users\DELL\Downloads'

    parser = argparse.ArgumentParser()
    parser.add_argument('--root', help='Root directory containing .svs files', default=def_root_dir)  # option that takes a value
    parser.add_argument('--clinical-data-json', help='Path to clinical data json downloaded from NCI GDC site', default=def_results_dir + '\\clinical.cohort.2024-04-04.json')
    parser.add_argument('--results', help='Results directory', default=def_results_dir)

    args = parser.parse_args()
    main(args)
