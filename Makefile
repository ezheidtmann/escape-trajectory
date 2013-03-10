
# Build tools and flags

CC			= build-tools/ccd-gcc
AR			= build-tools/run-ar
CFLAGS		= 

COMP		= @$(CC) $(CFLAGS) -Iframework -Iproblems -Iinclude -Iinclude/sundials $(TGT_CF) -o $@ -c $<
LINK		= @CC=g++ $(CC) $(CFLAGS) -o $@ $^ -Llib -lsundials_cvode -lsundials_nvecserial -lm
ARC			= @$(AR) -cr $@ $?

### 

.SUFFIXES:	.c .o

all: 		targets 

# subdirectories (like recursive make but better)

DIR := problems
include $(DIR)/Rules.mk
DIR := main
include $(DIR)/Rules.mk

# generic rules 

%.o:		%.c
	$(COMP)

%.oo:		%.cpp
	$(COMP)

%:			%.o
	$(LINK)

%:			%.oo
	$(LINK)

%:			*.a
	$(LINK)

escape: rthree
	cp rthree escape

rthree: $(OBJS) main/main-rthree.oo
	$(LINK)

twobody: $(OBJS) main/main-twobody.oo
	$(LINK)

curvedtwobody: $(OBJS) main/main-curved.oo
	$(LINK)

lyapunov-escape: $(OBJS) main/lyapunov-escape.oo
	$(LINK)

#$(TGTX): $(OBJS)
#	$(LINK) 

.PHONY:		targets
targets:	rthree escape twobody curvedtwobody $(TGTX)

CLEAN	:= $(CLEAN) $(TGTX)

.PHONY:		clean
clean:		
	rm -f $(CLEAN)

# vim:ts=4
