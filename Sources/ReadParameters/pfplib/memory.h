/* Copyright (C) 2003, 2004 Dave Bayer. Subject to the terms and conditions of the MIT License. */
#ifndef MEMORY_H
#define MEMORY_H
/*. MEMORYDEBUG */
#define MEMORYDEBUG1 1
#define MEMORYDEBUG 1

/*. memoryVerify */
#define memoryVerify memoryVerifyDebug

/*. malloc */
#define malloc(X) mallocDebug( (X), __FILE__, __LINE__, 0 )

/*. free */
#define free(X) freeDebug( (X), __FILE__, __LINE__ )

/*. mallocCount */
#define mallocCount mallocCountDebug

/* function prototypes */

void memoryVerifyDebug( void );
void *mallocDebug( size_t size, char *fileName, int lineNo, int ours );
void freeDebug( void *p, char *fileName, int lineNo );
void mallocCountDebug( void );

#endif
