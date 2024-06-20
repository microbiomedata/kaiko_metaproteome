# conda activate casanovo

FILES="../Kaiko_volume/Kaiko_input_files/mgf_large_unit_test/*.mgf"
for filepath in $FILES
do
    echo "Processing $filepath file..."
    arr=(${filepath//\// })
    mgf_filename=${arr[4]}
    
    echo 'casanovo_output\'"$mgf_filename"'.mztab'

    casanovo --mode=denovo --peak_path="$filepath" --output='casanovo_output\'"$mgf_filename"'.mztab'
done
