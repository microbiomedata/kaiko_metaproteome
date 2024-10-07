FROM camiloposso15/tensorflow2.12.0-py310

WORKDIR /Kaiko_metaproteome
ADD ./Kaiko_Deepnovo /Kaiko_metaproteome/Kaiko_Deepnovo
ADD ./Kaiko_main /Kaiko_metaproteome/Kaiko_main
ADD ./Casanovo /Kaiko_metaproteome/Casanovo
ADD ./Kaiko_dms /Kaiko_metaproteome/Kaiko_dms
ADD ./Kaiko_database_parse /Kaiko_metaproteome/Kaiko_database_parse
ADD ./Kaiko_unit_test /Kaiko_metaproteome/Kaiko_unit_test
ADD ./kaiko_defaults.yaml /Kaiko_metaproteome

RUN pip install --no-cache-dir --upgrade pip && \ 
    pip install --no-cache-dir indexed_gzip==1.8.5 llvmlite==0.39.1 biopython==1.81 numba==0.56.4 pyteomics sigopt==3.2.0 memory-profiler pyyaml pathlib s3path pyodbc openpyxl xlsxwriter torch torchvision torchaudio casanovo==4.2.1
