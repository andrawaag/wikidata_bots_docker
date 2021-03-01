FROM python:3.8-slim-buster
COPY . /app
WORKDIR /app
RUN apt-get update \
&& apt-get install gcc -y \
&& apt-get install git -y \
&& apt-get clean \
&& git clone https://github.com/sulab/scheduled-bots.git \
&& cd scheduled-bots && python3 setup.py install && pip install wikidataintegrator && pip install shexer && pip install --user -r requirements.txt

ENTRYPOINT ["python3", "scheduled-bots/scheduled_bots/geneprotein/GeneBot.py"]