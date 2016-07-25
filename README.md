# Source Code accompanying the paper "Unity and disunity in evolutionary sciences"

The paper "Unity and disunity in evolutionary sciencies" (List et al. 2016) introduces the concept of process-based analogies as a tool for interdisciplinary research. This repository contains python source code that can be used to replicate the examples on the use of similarity networks to detect borrowings and compound words. The output of the scripts are network-files in GML formats. In order to use these to create the images we provide in the paper, additional software is needed (we recommend Cytoscape, http://cytoscape.org).

## Requirements

* Python3
* lingpy (http://lingpy.org)
* networkx (http://networkx.github.io)

## Usage

There is one "library" script called `similarity_networks.py` and two scripts that create the networks:

* `Borrowing.py`, and 
* `Partial.py`

Just use them in the shell by calling:

```bash
$ python3 Borrowing.py
```

to create the file `r_person.gml`,

or 

```bash
$ python3 Partial.py
```

to create the file `r_face.gml`. All GML files are preceded by the prefix `r_`, indicating that this is part of the output produced by the scripts.
