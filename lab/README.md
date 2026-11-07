

O objetivo deste programa não é fornecer um solver completo, mas oferecer um ambiente controlado para estudar experimentalmente o comportamento de diferentes componentes do Multigrid.
Faça make help para obter instruções de uso.
Não se deve coletar os dados a partir de execuções isoladas, mas levantamento estatístico deve ser feito a partir do uso atento do código. Cabe ao usuário fazer esse estudo.
Teste várias configurações, mas também parâmetros diferentes do problema

leia problemas.md para ver limitações conhecidas e documentadas e extensões.md para ver extensões naturais a serem implementadas



- `lab.c`: menus e execução das combinações; `multigrid.c`: rotina FMG.
- `configuracao.*`: macros, configuração e nomes dos métodos; `resultados.h`: tipo dos resultados.
- `matriz.*`, `operadores.*`, `problema.*`, `residuo.*`, `erro.*`: malha e equação.
- `suavizacao.*`, `vcycle.*`, `direto.*`: algoritmos existentes, separados sem mudar fórmulas.
- `benchmark.*`: contadores e impressão do perfil da execução.

Para compilar: `make` (objetos `.o` + linkedição). Para os plots: `make dplot`; para o perfil: `make gprof`.
`make clean` remove objetos e executáveis.

Pesos padrao dos metodos: `OMEGA = 4.0/5.0` para Jacobi e `OMEGAGS = 1.0` para Gauss-Seidel.
Os dois valores estao em `configuracao.h` e sao aplicados automaticamente ao escolher as opcoes padrao nos blocos A/B/C e no modo esgotar.
No programa avulso `suavizador.c`, as duas macros estao definidas no proprio arquivo.
Os pesos continuam configuraveis quando o usuario opta por nao usar os padroes.
