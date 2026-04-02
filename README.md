![logo](nmrpf_logo_full.png)

[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/INRA/NMRProcFlow/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/INRA/NMRProcFlow)
[![NPFLOW](https://img.shields.io/badge/biotools-nmrprocflow-blue)](https://bio.tools/nmrprocflow)

## NMRProcFlow

An efficient GUI tool for spectra processing from 1D NMR metabolomics data

### Description

* The NMRProcFlow open source software provides an efficient GUI tool for spectra processing from 1D NMR metabolomics data, based on an interactive interface for the spectra visualization, that greatly helps spectra processing. 
* The 'NMR spectra viewer' is the central tool of NMRProcFlow and the core of the application. It allows the user to visually explore the spectra overlaid or stacked, to zoom on intensity scale, to group set of spectra by color based on their factor level.
* NMRProcFlow was built by involving NMR spectroscopists eager to have a quick and easy tool to use.
* the spectra processing includes: the calibration of the ppm scale, the base line correction (locally or  fully), the realignment of a set of ppm areas, and the binning (Intelligent, variable size or fixed width)


### Installation

Requirements:

* a recent Linux OS that support Docker (see https://www.docker.com/)


From GitHub, clone the repository

```
    $ git clone https://github.com/INRA/NMRProcFlow.git
```

Then `cd` to your clone path

```
    $ cd nmrproc
```

You need to have the docker image 'npflow' into your local docker library. For that, two ways are possible:

* You have to create the docker image:

		$ sh ./npflow build


* You can simply pull the docker image from docker hub:

		$ docker pull docker.io/nmrprocflow/npflow

Then, you should have something like below:

	$ docker image
```
REPOSITORY                   TAG                 IMAGE ID            CREATED             SIZE
nmrprocflow/nmrprocflow      latest              c7bf03be9a38        48 seconds ago      1.449 GB
```

Then edit **./etc/npflow.conf**, and pay attention to these lines :

```
# The URL root of the PROXY if applicable
PROXY_URL_ROOT=

# HTTP port number of the web application
HTTP_PORT=8080

# the root of the directory where data will be stored within a subdirectory and used by all modules
# The default path is /opt/data. If you choose to define a different location of the working space, 
# the corresponding directory will be created if it does not exist.
DATASETS=/opt/data

# Duration (in days) of validity of a session 
# before its destruction (counted from the last change)
PURGESESSIONS=2

# Max ZIP size (Mo)
MAXZIPSIZE=400

# NB CORES (0 means Auto)
CORES=0

# User connexion management
# 0 : no connexion management
# 1 : connexion management based on the /opt/data/conf/userlist file
#     Its structure is one user per line and each line following the format:
#        login;LastName;FirstName;Country;Institution;Email;Password
#     a minimal set of this 'userlist' file could be: npflow;;;;;;nppass
USRCONMGR=0
#
```


### Usage

#### Run the application

```
    $ sh ./npflow start
```

#### Stop the application

```
    $ sh ./npflow stop
```


#### View the status of the application

```
    $ sh ./npflow ps
```


#### check on http://<your_local_host>:<port>/npflow/


### Running on WSL2 (Windows Subsystem for Linux)

If you are running on WSL2 with Docker installed via **snap**, the default `/opt/data` path may be on a read-only filesystem. Follow these steps instead:

#### 1. Install Docker (snap)

If Docker is not yet installed:

```
    $ sudo snap install docker
    $ sudo snap start docker
```

Verify it works:

```
    $ sudo docker info
```

#### 2. Configure for WSL2

Edit `./etc/config.ini` and set `DEV=1` (development mode), then add a `DEV_DATADIR` line pointing to a writable location on your host:

```
    DEV=1
    HTTP_PORT=8081
    DATASETS=/opt/data
    DEV_DATADIR=/home/<your_username>/nmrprocflow_data
```

`DATASETS=/opt/data` is the path inside the container (do not change it). `DEV_DATADIR` is the host-side directory used as the Docker volume mount source.

Create the data directory:

```
    $ mkdir -p /home/<your_username>/nmrprocflow_data
```

#### 3. Build and run

```
    $ sudo sh ./npflow build
    $ sudo sh ./npflow start
```

In dev mode (DEV=1), the application runs on the port specified in `config.ini` (default 8081). Open your browser at:

```
    http://localhost:8081/npflow/
```

Dev mode mounts `./nmrspec/` and `./nmrviewer/www/` into the container, so code changes on your host are reflected without rebuilding.


### Stopping and cleaning up

#### Stop the application

```
    $ sh ./npflow stop
```

This stops and removes the running container.

#### Remove the Docker image

```
    $ sudo docker rmi nmrprocflow/nmrprocflow:latest
```

#### Remove all data

```
    $ rm -rf /home/<your_username>/nmrprocflow_data
```

Or if you used the default path:

```
    $ sudo rm -rf /opt/data
```

#### Full cleanup (stop container, remove image, prune)

```
    $ sh ./npflow stop
    $ sudo docker rmi nmrprocflow/nmrprocflow:latest
    $ sudo docker system prune -f
```

This removes the container, image, and any dangling Docker resources.

#### Removing the Docker bridge network (WSL2)

On some WSL2 networks, the `docker0` bridge interface can interfere with routing (e.g. SSH connections failing with "No route to host"). This happens when Docker's default bridge subnet overlaps with your destination network.

**Important:** Docker needs the `docker0` bridge to build images and run containers. Only remove it after you are done using Docker.

Recommended workflow:

```
    # 1. Start Docker, build, and run
    $ sudo snap restart docker
    $ sudo sh ./npflow build
    $ sudo sh ./npflow start

    # 2. When done, stop everything and remove the bridge
    $ sudo sh ./npflow stop
    $ sudo snap stop docker
    $ sudo ip link set docker0 down
    $ sudo ip link delete docker0
```

To prevent Docker from creating the bridge at all (note: this breaks `docker build`):

```
    $ sudo mkdir -p /etc/docker
    $ echo '{"bridge": "none"}' | sudo tee /etc/docker/daemon.json
```

If you use this option, remove `/etc/docker/daemon.json` before building again.


See more information in http://nmrprocflow.org/c_download


### Project structure

```
nmrspec/
├── server.R                  # Shiny server entry point
├── global.R                  # Global setup, sources RnmrTools.R
├── exec/
│   ├── RnmrTools.R           # Main loader — sources all modules below
│   ├── lib/
│   │   ├── utils.R           # String helpers, counters, INI/LOG I/O, default_noise_range()
│   │   ├── stack.R           # Undo history (push/pop/clean)
│   │   ├── metadata.R        # Sample metadata generation from raw spectra archives
│   │   ├── baseline.R        # Baseline correction (Whittaker, airPLS, global, local, q-NMR)
│   │   ├── alignment.R       # Spectral alignment (CluPA, PTW, least-squares, shift)
│   │   ├── normalization.R   # Calibration and normalization (CSN, PQN)
│   │   ├── processing.R      # Denoising (Savitzky-Golay), zeroing, smoothing
│   │   ├── bucketing.R       # Spectral binning algorithms and macro-command validation
│   │   ├── macro.R           # Macro-command file replay (RProcCMD1D)
│   │   └── export.R          # Bucket/SNR dataset generation and spectra data export
│   ├── libspec/              # C++ (Rcpp) low-level spectral functions
│   ├── Rcorr1D               # Preprocessing entry point
│   ├── Ralign1D              # Alignment entry point
│   ├── Rbuc1D                # Bucketing entry point
│   └── Rnmr1D                # Main processing entry point
├── Rsrc/
│   ├── Proc1.R – Proc4.R    # Shiny server modules (upload, processing, bucketing, export)
│   ├── ui_procparams.R       # Processing parameter UI
│   └── utils.R               # Shiny utility functions
└── conf/                     # Apache, Shiny Server, and app configuration
```

`RnmrTools.R` is a thin loader that imports libraries and sources the 10 modules in `exec/lib/`. All callers (`global.R`, `Rnmr1D`, `Rcorr1D`, `Ralign1D`, `Rbuc1D`) continue to `source("exec/RnmrTools.R")` — no caller changes needed.

---

### Main contributors:

* Daniel Jacob, Catherine Deborde, Marie Lefebvre, Michaël Maucourt
* Special thanks to Alain Girard (INRA Bordeaux) for designing the logo.

### Publication:

Jacob, D., Deborde, C., Lefebvre, M., Maucourt, M. and Moing, A. (2017) NMRProcFlow: A graphical and interactive tool dedicated to 1D spectra processing for NMR-based metabolomics, Metabolomics 13:36. [doi:10.1007/s11306-017-1178-y](http://link.springer.com/article/10.1007%2Fs11306-017-1178-y)

### Funded by:

* INRAE UMR 1332 BFP, Bordeaux Metabolomics Facility
* the ANR-11-INBS-0010 grant (MetaboHUB)

### License

Copyright (C) 2016-2021  Daniel Jacob - INRAE 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
