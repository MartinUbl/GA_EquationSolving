#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "stack.h"

c_stack* stck_create(int size)
{
    c_stack *stck = new c_stack();

    stck->size = size;
    stck->curr = STCK_INVALID;
    stck->elements = (void**)malloc(sizeof(void*)*size);
    memset(stck->elements, 0, sizeof(void*)*size);

    return stck;
}

void stck_destroy(c_stack *stck)
{
    /* not popped elements will be deleted */
    for (; stck->curr >= 0; stck->curr--)
    {
        if (stck->elements[stck->curr] != NULL)
            free(stck->elements[stck->curr]);
    }

    free(stck->elements);
    delete stck;
}

void* stck_pop(c_stack *stck)
{
    /* STCK_INVALID is -1, so it's set automatically when the stack is emptied */
    if (stck->curr == STCK_INVALID)
        return NULL;

    return stck->elements[stck->curr--];
}

void* stck_peek(c_stack *stck)
{
    return stck_get(stck, stck->curr);
}

void* stck_get(c_stack *stck, int pos)
{
    /* STCK_INVALID is -1, so it's set automatically when the stack is emptied */
    if (stck->curr == STCK_INVALID)
        return NULL;

    return stck->elements[pos];
}

void stck_push(c_stack *stck, void *el)
{
    if (stck->curr == stck->size - 1)
        return;

    stck->elements[++stck->curr] = el;
}
