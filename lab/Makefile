CC = gcc
CFLAGS = -Wall -Wextra -O2
EXTRA_CFLAGS =
LDLIBS = -lm

include config.mk

TARGET = lab
MODULOS = lab configuracao benchmark matriz problema operadores residuo suavizacao \
          vcycle direto erro multigrid
OBJS = $(addsuffix .o,$(MODULOS))
PLOT_OBJS = $(addsuffix .plot.o,$(MODULOS))
PROF_OBJS = $(addsuffix .prof.o,$(MODULOS))
VIDEO_OBJS = $(addsuffix .video.o,$(MODULOS))
VIDEO_N ?= 64

ifeq ($(PRINT), ON)
CFLAGS += -DPRINT
endif

.PHONY: all clean config unprint print profile dplot plot gprof help video

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $^ -o $@ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -c $< -o $@

# Objetos de plot separados para manter a compilação normal intacta.
dplot: labplot
labplot: $(PLOT_OBJS)
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $^ -o $@ $(LDLIBS)
%.plot.o: %.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -DPLOT -c $< -o $@

# Objetos de perfil separados; -pg também é usado na linkedição.
gprof: labprof
labprof: $(PROF_OBJS)
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -pg $^ -o $@ $(LDLIBS)
%.prof.o: %.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -pg -c $< -o $@

plot:
	@bash plot.sh

clean:
	rm -f $(TARGET) labplot labprof profile *.o *.txt *.out *.zip *.vtk *.dat

config:
	@cat config.mk

print:
	@sed -i 's/^PRINT = .*/PRINT = ON/' config.mk
	@rm -f $(TARGET) *.o
	@$(MAKE) config

unprint:
	@sed -i 's/^PRINT = .*/PRINT = OFF/' config.mk
	@rm -f $(TARGET) *.o
	@$(MAKE) config

profile:
	@bash profile.sh

# Compilacao independente, sem contaminar os objetos normais.
labvideo: $(VIDEO_OBJS)
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -DVIDEO $^ -o $@ $(LDLIBS)
%.video.o: %.c
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) -DVIDEO -c $< -o $@

video:
	@bash video.sh "$(VIDEO_N)"

help:
	@printf '%s\n' \
	  'LABORATORIO MULTIGRID - alvos do Makefile' \
	  '  make                 compila o executavel normal lab' \
	  '  make help            mostra esta ajuda' \
	  '  make video           cria video.mp4 3D com malha N=64 (padrao)' \
	  '  make video VIDEO_N=128   escolhe outra malha: 16,32,64,128,...' \
	  '                       pergunta parametros A/B/C na execucao' \
	  '                       usa gnuplot, ffmpeg e ffplay (-loop 0)' \
	  '                       remove labvideo, objetos e quadros ao terminar' \
	  '  VIDEO_NO_PLAY=1 make video    gera MP4 sem abrir ffplay' \
	  '  make dplot          compila labplot com -DPLOT' \
	  '  make plot           roda plot.sh (graficos estaticos)' \
	  '  make gprof          compila labprof com -pg' \
	  '  make profile        roda profile.sh' \
	  '  make print          ativa PRINT em config.mk e recompila ao usar make' \
	  '  make unprint        desativa PRINT em config.mk e recompila ao usar make' \
	  '  make config         mostra config.mk' \
	  '  make clean          remove artefatos conforme regra original' \
	  '' \
	  'Execucao comum: ./lab N [modo]; N deve ser 16*2^k.' \
	  'Com N sozinho, o programa pergunta se deseja escolher configuracao.' \
	  'Modo 0: esgota 50 combinacoes A/B/C; modo 1: escolha interativa.' \
	  'ATENCAO: make clean usa a regra original, que apaga *.txt e *.zip.'
