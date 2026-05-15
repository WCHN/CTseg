FROM ubuntu:22.04

LABEL org.opencontainers.image.authors="SPM <fil.spm@ucl.ac.uk>"
LABEL org.opencontainers.image.source="https://github.com/WCHN/CTseg"

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
     unzip xorg wget \
 && apt-get clean \
 && rm -rf \
     /tmp/hsperfdata* \
     /var/*/apt/*/partial \
     /var/lib/apt/lists/* \
     /var/log/apt/term*

# Install MATLAB MCR in /opt/mcr/
ENV MATLAB_VERSION=R2026a
ENV MCR_DIR=R2026a
RUN mkdir /opt/mcr_install \
 && mkdir /opt/mcr \
 && wget --progress=bar:force -P /opt/mcr_install https://ssd.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_${MATLAB_VERSION}_glnxa64.zip \
 && unzip -q /opt/mcr_install/MATLAB_Runtime_${MATLAB_VERSION}_glnxa64.zip -d /opt/mcr_install \
 && /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent \
 && rm -rf /opt/mcr_install /tmp/*

# Install SPM Standalone in /opt/spm12/
ENV SPM_REVISION=r7771
ENV LD_LIBRARY_PATH=/opt/mcr/${MCR_DIR}/runtime/glnxa64:/opt/mcr/${MCR_DIR}/bin/glnxa64:/opt/mcr/${MCR_DIR}/sys/os/glnxa64:/opt/mcr/${MCR_DIR}/sys/opengl/lib/glnxa64:/opt/mcr/${MCR_DIR}/extern/bin/glnxa64
ENV MCR_INHIBIT_CTF_LOCK=1
ENV SPM_HTML_BROWSER=0
# Running SPM once with "function exit" tests the succesfull installation *and*
# extracts the ctf archive which is necessary if singularity is going to be
# used later on, because singularity containers are read-only.
# Also, set +x on the entrypoint for non-root container invocations
ADD https://github.com/WCHN/CTseg/releases/download/v1.0.2/spm12_${SPM_REVISION}_BI_Linux_${MATLAB_VERSION}.zip /opt
RUN unzip -q /opt/spm12_${SPM_REVISION}_BI_Linux_${MATLAB_VERSION}.zip -d /opt \
 && rm -f /opt/spm12_${SPM_REVISION}_BI_Linux_${MATLAB_VERSION}.zip \
 && /opt/spm12/spm12 function exit \
 && chmod +x /opt/spm12/spm12

# Download default atlas to models/ directory on first use
RUN /opt/spm12/spm12 eval "try,spm_CTseg(1);end"

# Configure entry
ENTRYPOINT ["/opt/spm12/spm12"]
CMD ["--help"]