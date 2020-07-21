FROM python:3.7
RUN mkdir /sandbox /input_files /resources /output_files /calculate_common_coverage
ADD calculate_common_coverage /calculate_common_coverage
ADD sambamba_0.7.1 /resources
RUN ls /resources
WORKDIR /calculate_common_coverage
RUN python3 /calculate_common_coverage/calculate_common_coverage.py -h
ENTRYPOINT [""]

