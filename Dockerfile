FROM tensorflow2.12.0-py310

WORKDIR /Kaiko_pipeline
ADD ./Kaiko_denovo /Kaiko_pipeline/Kaiko_denovo
ADD ./Kaiko_2.py /Kaiko_pipeline
ADD ./Kaiko_3.py /Kaiko_pipeline
ADD ./Kaiko_4.py /Kaiko_pipeline
ADD ./Kaiko_pipeline_main.py /Kaiko_pipeline
ADD ./kaiko_unit_test.py /Kaiko_pipeline
ADD ./unit_test_util.py /Kaiko_pipeline
ADD ./kaiko_defaults.yaml /Kaiko_pipeline

RUN pip install --no-cache-dir --upgrade pip && \ 
    pip install --no-cache-dir llvmlite==0.39.1 biopython==1.81 numba==0.56.4 pyteomics sigopt==3.2.0 memory-profiler pyyaml pathlib s3path
