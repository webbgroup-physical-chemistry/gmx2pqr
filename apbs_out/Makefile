NAME := read_apbs_rxn_field
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

CXX := g++ 
CXXFLAGS := -O3 -Wall 
INCLUDE := 
LDFLAGS :=
LDLIBS :=

$(NAME): $(OBJS)
	$(CXX) $(LDFLAGS) $(LDLIBS) -o $@ $^
$.o: %.c
	$(CXX) $(CXXFLAGS) -c $<
clean:
	rm -f $(NAME) $(OBJS)
print-%:
	@echo '$* =  $($*)'

.PHONY: clean
