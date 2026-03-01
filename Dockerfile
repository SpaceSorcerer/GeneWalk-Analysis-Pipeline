FROM python:3.11-slim

WORKDIR /app

# Install web visualizer dependencies (no GeneWalk needed — web apps only
# visualize pre-computed results or run GSEA/ORA which are pure-Python)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

HEALTHCHECK CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health')" || exit 1

# Default: GeneWalk results visualizer.
# To run the comparison pipeline instead:
#   docker run -p 8501:8501 genewalk-app comparison_app.py
ENTRYPOINT ["streamlit", "run", "--server.port=8501", "--server.address=0.0.0.0"]
CMD ["app.py"]
