# DF2D_Altos
DART-FISH analysis pipeline for internal use at Altos

To use the pipeline:
1. Install the main conda environment: `conda env create -f conda_environment_image.yml`
2. Activate conda environment: `conda activate image`
3. Modify parameter file for your experiment.
    Important parameters are: `reg_rounds`, `rounds_dir`, `rounds_regex3d` for file structures; `ref_reg_cycle` for cycles as reference for registration; `codebook_path`, `dc_rounds` for decoding; `pretrained_model` for segementation
4. Run the pipeline in dockyard using: `python main.py params_2dreg.yaml`
5. Or to run the pipeline step by step:
    ```
    python file_adapter.py params_2dreg.yaml
    python rawMaxProjection.py params_2dreg.yaml
    python AlignerDriver_2D.py params_2dreg.yaml
    python StitchDriver_2D.py params_2dreg.yaml
    python findDecodingParams.py params_2dreg.yaml
    python decoding_driver.py params_2dreg.yaml
    # if not enough memory, use docker job --type r6i.32xlarge --cpu 15 --memory 500 python_path decoding_driver.py params_2dreg.yaml
    python CombineFOVs.py params_2dreg.yaml
    python segmentation_driver.py params_2dreg.yaml
    python QC_plots.py params_2dreg.yaml
    ``` 
