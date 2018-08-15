conda alt install

For users without a home directory on the genome center clusters, this document describes
alternate installation of conda that doesn't require a home folder. You won't be able to
build environments, but you can use Kirk's without having to modifying any additional code.

1. Download miniconda

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2. Run miniconda installation script, but change home location to a folder you can access,
e.g., ```/bruno/comai/yourUserName/miniconda3```

```
bash Miniconda3-latest-Linux-x86_64.sh
# review license
# accept license
# change home location to something like /bruno/comai/yourUserName/miniconda3
# no to placing it in your path
```

3. Append conda location to your path, e.g.,

```
export PATH='/bruno/comai/yourUserName/miniconda3/bin:$PATH'
```

4. Check that this worked by using the command ```which conda```. This should return

```/bruno/comai/yourUserName/miniconda3/bin/conda```
