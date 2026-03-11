FROM python:3.11-slim

WORKDIR /app

# Install web visualizer dependencies (no GeneWalk needed -- web app only
# visualizes pre-computed results unless genewalk is installed separately)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

HEALTHCHECK CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

ENTRYPOINT ["streamlit", "run", "--server.port=8501", "--server.address=0.0.0.0"]
CMD ["app.py"]
