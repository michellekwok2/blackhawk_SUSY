TARGET    = fM_DM2.x
OBJECT_01 = fM_DM2.o
SOURCE_01 = fM_DM2.c

$(TARGET): $(OBJECT_01)
	gcc -g -o $(TARGET) $(OBJECT_01) -lm -std=c99

$(OBJECT_01) : $(SOURCE_01)
	gcc -g -c $(SOURCE_01) -o $(OBJECT_01) -lm -std=c99

all : $(TARGET)

clean :
	-rm $(TARGET) $(OBJS)
