FROM python:3.12-slim

# gcc + make: o app compila o executável `lab` em tempo de execução
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc make libc6-dev \
 && rm -rf /var/lib/apt/lists/*

ENV PYTHONUNBUFFERED=1 MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl
WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Pré-compila os executáveis para EXP = 1, 2, 3 (evita compilar na 1ª requisição)
RUN python -c "import app; [app._garante_executaveis(e) for e in (1, 2, 3)]"

RUN useradd -m appuser && chown -R appuser /app
USER appuser

# Render injeta PORT (padrão 10000); o servidor precisa escutar em 0.0.0.0
CMD ["sh", "-c", "gunicorn app:app --bind 0.0.0.0:${PORT:-10000} --timeout 120 --threads 1"]
