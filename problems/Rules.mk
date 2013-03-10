
# Standard stuff

sp			:= $(sp).x
dirstack_$(sp) := $(d)
d			:= $(DIR)

# add subdirectories here

OBJS_$(d)	:= $(d)/RestrictedThree.oo $(d)/AbstractProblem.oo $(d)/SingleBody.oo $(d)/TwoBody.oo $(d)/CurvedTwoBody.oo

OBJS		:= $(OBJS) $(OBJS_$(d))
DEPS_$(d)	:= $(OBJS_$(d):%=%.d)
CLEAN		:= $(CLEAN) $(OBJS_$(d)) $(DEPS_$(d)) 

# this inherits from other stuff to do the work
$(OBJS_$(d)): TGT_CF := -I$(d) -Ilibserial

# include dependency files
-include 	$(DEPS_$(d))

# pop off the stack

d 			:= $(dirstack_$(sp))
sp			:= $(basename $(sp))

# vim:ts=4
