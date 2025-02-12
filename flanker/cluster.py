import os
import glob
import tempfile
import subprocess
import collections
import pandas as pd
import networkx as nx

"""
functions to cluster output sequences from flanker
build mash sketch and pairwise_mash_distances are functions written by Ryan Wick in his assembly de-replicator repo
https://github.com/rrwick/Assembly-Dereplicator. The other functions are also adapted from functions written by Ryan
for assembly de-replicator.
"""

def find_all_assemblies(output_dir):
    all_assemblies = []
    missing_fastas = []

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        print(f"Warning: Output directory {output_dir} does not exist.")
        return all_assemblies

    # Only search inside the specified output directory
    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)

        # Print the file path to check what it's finding
        print(f"Checking file: {file_path}")

        if os.path.isfile(file_path) and os.path.getsize(file_path) > 100 and filename.endswith('flank.fasta'):
            try:
                with open(file_path, 'r') as f:
                    pass  # Just attempting to open it
                all_assemblies.append(file_path)
            except Exception:
                missing_fastas.append(file_path)

    # Log missing files
    if missing_fastas:
        with open("missing_fastas.txt", "a") as log_file:
            log_file.write("\n".join(missing_fastas) + "\n")

    print(f'Found {len(all_assemblies):,} files in {output_dir}')
    print(f'Logged {len(missing_fastas)} missing files to missing_fastas.txt')

    return all_assemblies

def build_mash_sketch(assemblies, threads, temp_dir, sketch_size,kmer_length):
    mash_command = ['mash', 'sketch', '-p', str(threads), '-o', temp_dir + '/mash',
                        '-s', str(sketch_size), '-k', str(kmer_length)] + assemblies


    subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return temp_dir + '/mash.msh'

def pairwise_mash_distances(mash_sketch, threads):
    mash_command = ['mash', 'dist', '-p', str(threads), mash_sketch, mash_sketch]

    mash_out = subprocess.run(mash_command, stdout=subprocess.PIPE).stdout.decode()
    return mash_out.splitlines()



def create_graph_from_distances(pairwise_distances, threshold):

    matrix=[]
    assemblies = set()
    graph = collections.defaultdict(set)
    all_connections = collections.defaultdict(set)
    for line in pairwise_distances:
        parts = line.split('\t')
        assembly_1 = parts[0]
        assembly_2 = parts[1]
        distance = float(parts[2])

        matrix.append([assembly_1,assembly_2,distance])



    df=pd.DataFrame(matrix,columns=['assembly_1','assembly_2','distance'])
    df=df[df.distance <= threshold]
    G=nx.from_pandas_edgelist(df,'assembly_1','assembly_2','distance')
    l=list(nx.connected_components(G))
    L=[dict.fromkeys(y,x) for x, y in enumerate(l)]

    d={k: v for d in L for k, v in d.items()}
    df2=df['assembly_1'].unique()
    df2=pd.DataFrame(df2,columns=['assembly_1'])

    df2['cluster']=df2.assembly_1.map(d)

    return(df2)

def flank_scrub():

    filelist=glob.glob(str(str(os.getcwd()) + '/' + str("*flank.fasta")))

    for filename in filelist:
        os.remove(filename)


#here we build clusters using mash distances
def define_clusters(gene, window, threads, threshold, outfile, kmer_length, sketch_size):
    with tempfile.TemporaryDirectory() as temp_dir:
        all_assemblies = find_all_assemblies(outfile)  # Pass the correct directory
        mash_sketch = build_mash_sketch(all_assemblies, threads, temp_dir, sketch_size, kmer_length)

        pairwise_distances = pairwise_mash_distances(mash_sketch, threads)
        clusters = create_graph_from_distances(pairwise_distances, float(threshold))

        clusters.to_csv(f"{outfile}_{gene.strip()}_{window}", index=False)
