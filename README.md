# Laboratório Multigrid na web (Flask)

    pip install -r requirements.txt
    python3 app.py            # http://localhost:5000

Requer `gcc` e `make`. Os fontes C ficam em `lab/` (ou aponte com `LAB_DIR=/caminho/lab`).
O app não altera seus fontes nem seu `config.mk`: compila executáveis próprios em `_cache/EXP<n>/`
(EXP é macro de compilação; K, ω, ν1, ν2, tolerâncias, A, B e C são passados ao `lab` por stdin)
e roda cada requisição numa pasta temporária.

- **Executar**: roda `./labplot N 1` com a configuração dos sliders. O programa resolve N e 2N
  (para estimar a ordem p); o gráfico (analítica, numérica, erro) é o da malha 2N.
  Também mostra o histórico do resíduo, o perfil de tempo e a saída completa.
- **Comparar as 50 combinações**: equivalente ao `./lab N 0` (modo "esgotar"), com tabela
  ordenada por tempo e mapas de calor (B × C) de tempo, ciclos V, ordem p e operações na malha grossa.
  Cada método usa seu peso padrão; o restante segue os sliders.
- Limites: N ≤ 512 (execução simples), N ≤ 256 (comparação); ver `KMAX` e `KMAX_COMPARAR`.
- Tempos de execuções isoladas são ruidosos; como diz o README do lab, estatística exige repetições.
