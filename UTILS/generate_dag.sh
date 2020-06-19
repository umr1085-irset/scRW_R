#!/bin/bash
function yaml() {
    hashdot=$(gem list hash_dot);
    if ! [ "$hashdot" != "" ]; then sudo gem install "hash_dot" ; fi
    if [ -f $1 ];then
        cmd=" Hash.use_dot_syntax = true; hash = YAML.load(File.read('$1'));";
        if [ "$2" != "" ] ;then 
            cmd="$cmd puts hash.$2;"
        else
            cmd="$cmd puts hash;"
        fi
        ruby -r yaml -r hash_dot <<< $cmd;
    fi
}
outdir=$(yaml config.yaml outdir) # get outdir value from config file
outdirdag=$outdir/DAG
mkdir -p $outdirdag # create outdir folder if necessary
snakemake --dag | dot -Tsvg > $outdirdag/dag.svg # save dag to svg
snakemake --dag | dot -Tpng > $outdirdag/dag.png # save dag to png