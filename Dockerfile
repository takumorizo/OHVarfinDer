FROM centos:6

RUN set -x && \
    yum -y install gcc  && \
    yum -y install gcc-c++ && \
    yum -y install git     && \
    yum -y install ncurses-devel && \
    yum -y install zlib-devel    && \
    cd /opt/ && \
    git clone https://github.com/takumorizo/OHVarfinDer.git && \
    cd ./OHVarfinDer && \
    make dependencies && \
    make

ENTRYPOINT ["/opt/OHVarfinDer/bin/ohvarfinder"]
