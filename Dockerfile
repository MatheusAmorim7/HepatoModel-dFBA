# Hepatic lobule dFBA: reproducible environment (Ubuntu 24.04) with JupyterLab.
#   docker build -t lobule-dfba .
#   docker run --rm -it -p 8888:8888 lobule-dfba          # open the printed http://127.0.0.1:8888/lab?token=... link
#   docker run --rm -it lobule-dfba make quick             # headless smoke test (about 5-10 min)
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      g++ make git coinor-libclp-dev coinor-libcoinutils-dev libsbml5-dev \
      liblapack-dev libxml2-dev libbz2-dev zlib1g-dev python3 python3-venv python3-pip ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir numpy scipy pandas cobra matplotlib jupyterlab
ENV PATH="/opt/venv/bin:${PATH}"
WORKDIR /opt/hepatomodel
COPY . .
WORKDIR /opt/hepatomodel/sample_projects_intracellular/fba/hepatic_lobule_dfba
RUN make -j"$(nproc)" PY=python3
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--allow-root", "--notebook-dir=/opt/hepatomodel/sample_projects_intracellular/fba/hepatic_lobule_dfba/notebooks"]
