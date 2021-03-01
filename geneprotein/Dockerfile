FROM python:3.8-slim-buster
COPY . /app
WORKDIR /app
RUN apt-get update \
&& apt-get install gcc -y \
&& apt-get clean && pip install --user -r requirements.txt
RUN python3 setup.py install
CMD ["python3", "GeneBot.py"]