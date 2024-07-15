FROM mrtrix3/mrtrix3
FROM freesurfer/freesurfer:7.1.1
COPY license /usr/local/freesurfer/.license

ENV FSLDIR          "/usr/local/fsl"
ENV DEBIAN_FRONTEND "noninteractive"
ENV LANG            "en_GB.UTF-8"

RUN curl https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz \
      | tar -xz -C /usr/local && \
       /usr/local/fsl/etc/fslconf/fslpython_install.sh -f /usr/local/fsl

# Configure environment
ENV FSLDIR=/usr/local/fsl
ENV FSL_DIR="${FSLDIR}" \
    FSLOUTPUTTYPE=NIFTI_GZ \
    PATH=${FSLDIR}/bin:$PATH \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=${FSLDIR} \
    LD_LIBRARY_PATH=${FSLDIR}/lib:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    FSLOUTPUTTYPE=NIFTI_GZ

# Install FreeSurfer
RUN wget -qO- ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu20_amd64-7.4.1.tar.gz | tar -xz -C /opt \
&& rm -rf /opt/freesurfer/subjects \
&& rm -rf /opt/freesurfer/trctrain


# install gradient_unwarp.py (v1.2.0 with python 3 compatibility)
WORKDIR /tmp
RUN wget -q https://github.com/Washington-University/gradunwarp/archive/v1.2.0.zip && \
  unzip v1.2.0.zip && \
  cd gradunwarp-1.2.0 && \
  python setup.py install && \
  rm -rf gradunwarp-1.2.0 v1.2.0.zip

# Install MCR 2017b
ENV MATLABCMD="/opt/matlabmcr-2017b/v93/toolbox/matlab" \
    MATLAB_COMPILER_RUNTIME="/opt/matlabmcr-2017b/v93" \
    LD_LIBRARY_PATH="/opt/matlabmcr-2017b/v93/runtime/glnxa64:/opt/matlabmcr-2017b/v93/bin/glnxa64:/opt/matlabmcr-2017b/v93/sys/os/glnxa64:$LD_LIBRARY_PATH"

# overwrite matlab mcr shared object
RUN rm /opt/matlabmcr-2017b/v93/sys/os/glnxa64/libstdc++.so.6 && \
    ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /opt/matlabmcr-2017b/v93/sys/os/glnxa64/libstdc++.so.6

ENTRYPOINT [ "sh", "-c", ". /usr/local/fsl/etc/fslconf/fsl.sh && /bin/bash" ]