
# Standard stuff

sp			:= $(sp).x
dirstack_$(sp) := $(d)
d			:= $(DIR)

# add subdirectories here

OBJS_$(d)	:= $(d)/main-twobody.oo $(d)/main-rthree.oo $(d)/main-curved.oo $(d)/lyapunov-escape.oo

# disable adding these to the big list because we handle them directly here.
#OBJS		:= $(OBJS) $(OBJS_$(d))
DEPS_$(d)	:= $(OBJS_$(d):%=%.d)
CLEAN		:= $(CLEAN) $(OBJS_$(d)) $(DEPS_$(d)) 

$(d)/main-twobody.oo: $(d)/main.cpp
	$(COMP) -DPROBLEM_TYPE=TwoBody

$(d)/main-rthree.oo: $(d)/main.cpp
	$(COMP) -DPROBLEM_TYPE=RestrictedThree

$(d)/main-curved.oo: $(d)/main.cpp
	$(COMP) -DPROBLEM_TYPE=CurvedTwoBody

$(d)/lyapunov-escape.oo: $(d)/lyapunov-escape.cpp
	$(COMP)

# this inherits from other stuff to do the work
#$(OBJS_$(d)): TGT_CF := -I$(d) -Ilibserial

# include dependency files
-include 	$(DEPS_$(d))

# pop off the stack

d 			:= $(dirstack_$(sp))
sp			:= $(basename $(sp))

# vim:ts=4
