/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Only edit the mpe_thread.h.in version of this file:
 *  src/util/thread/mpe_thread.h.  Generated from mpe_thread.h.in by configure.
 */

#if !defined(MPE_THREAD_H_INCLUDED)
#define MPE_THREAD_H_INCLUDED

#if !defined(TRUE)
#define TRUE 1
#endif
#if !defined(FALSE)
#define FALSE 0
#endif


/*
 * Implementation specific type definitions
 */
/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <errno.h>
#include <pthread.h>

typedef pthread_mutex_t MPE_Thread_mutex_t;
typedef pthread_cond_t  MPE_Thread_cond_t;
typedef pthread_t       MPE_Thread_id_t;
typedef pthread_key_t   MPE_Thread_tls_t;

#define MPE_THREAD_TLS_T_NULL 0


/*
 * Threads
 */

typedef void (* MPE_Thread_func_t)(void * data);

/*@
  MPE_Thread_create - create a new thread

  Input Parameters:
+ func - function to run in new thread
- data - data to be passed to thread function

  Output Parameters:
+ id - identifier for the new thread
- err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred

  Notes:
  The thread is created in a detach state, meaning that is may not be waited upon.  If another thread needs to wait for this
  thread to complete, the threads must provide their own synchronization mechanism.
@*/
void MPE_Thread_create(MPE_Thread_func_t func, void * data, MPE_Thread_id_t * id, int * err);

/*@
  MPE_Thread_exit - exit from the current thread
@*/
void MPE_Thread_exit(void);

/*@
  MPE_Thread_self - get the identifier of the current thread

  Output Parameter:
. id - identifier of current thread
@*/
void MPE_Thread_self(MPE_Thread_id_t * id);

/*@
  MPE_Thread_same - compare two threads identifiers to see if refer to the same thread

  Input Parameters:
+ id1 - first identifier
- id2 - second identifier
  
  Output Parameter:
. same - TRUE if the two threads identifiers refer to the same thread; FALSE otherwise
@*/
void MPE_Thread_same(MPE_Thread_id_t * id1, MPE_Thread_id_t * id2, int * same);

/*@
  MPE_Thread_yield - voluntarily relinquish the CPU, giving other threads an opportunity to run
@*/
void MPE_Thread_yield(void);


/*
 *    Mutexes
 */

/*@
  MPE_Thread_mutex_create - create a new mutex
  
  Output Parameters:
+ mutex - mutex
- err - error code (non-zero indicates an error has occurred)
@*/
void MPE_Thread_mutex_create(MPE_Thread_mutex_t * mutex, int * err);

/*@
  MPE_Thread_mutex_destroy - destroy an existing mutex
  
  Input Parameter:
. mutex - mutex

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_mutex_destroy(MPE_Thread_mutex_t * mutex, int * err);

/*@
  MPE_Thread_lock - acquire a mutex
  
  Input Parameter:
. mutex - mutex

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_mutex_lock(MPE_Thread_mutex_t * mutex, int * err);

/*@
  MPE_Thread_unlock - release a mutex
  
  Input Parameter:
. mutex - mutex

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_mutex_unlock(MPE_Thread_mutex_t * mutex, int * err);

/*@
  MPE_Thread_mutex_trylock - try to acquire a mutex, but return even if unsuccessful
  
  Input Parameter:
. mutex - mutex

  Output Parameters:
+ flag - flag
- err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_mutex_trylock(MPE_Thread_mutex_t * mutex, int * flag, int * err);


/*
 * Condition Variables
 */

/*@
  MPE_Thread_cond_create - create a new condition variable
  
  Output Parameters:
+ cond - condition variable
- err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_cond_create(MPE_Thread_cond_t * cond, int * err);

/*@
  MPE_Thread_cond_destroy - destroy an existinga condition variable
  
  Input Parameter:
. cond - condition variable

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_cond_destroy(MPE_Thread_cond_t * cond, int * err);

/*@
  MPE_Thread_cond_wait - wait (block) on a condition variable
  
  Input Parameters:
+ cond - condition variable
- mutex - mutex

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred

  Notes:
  This function may return even though another thread has not requested that a thread be released.  Therefore, the calling
  program must wrap the function in a while loop that verifies program state has changed in a way that warrants letting the
  thread proceed.
@*/
void MPE_Thread_cond_wait(MPE_Thread_cond_t * cond, MPE_Thread_mutex_t * mutex, int * err);

/*@
  MPE_Thread_cond_broadcast - release all threads currently waiting on a condition variable
  
  Input Parameter:
. cond - condition variable

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_cond_broadcast(MPE_Thread_cond_t * cond, int * err);

/*@
  MPE_Thread_cond_signal - release one thread currently waitng on a condition variable
  
  Input Parameter:
. cond - condition variable

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_cond_signal(MPE_Thread_cond_t * cond, int * err);


/*
 * Thread Local Storage
 */
typedef void (*MPE_Thread_tls_exit_func_t)(void * value);


/*@
  MPE_Thread_tls_create - create a thread local storage space

  Input Parameter:
. exit_func - function to be called when the thread exists; may be NULL is a callback is not desired
  
  Output Parameters:
+ tls - new thread local storage space
- err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_tls_create(MPE_Thread_tls_exit_func_t exit_func, MPE_Thread_tls_t * tls, int * err);

/*@
  MPE_Thread_tls_destroy - destroy a thread local storage space
  
  Input Parameter:
. tls - thread local storage space to be destroyed

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
  
  Notes:
  The destroy function associated with the thread local storage will not called after the space has been destroyed.
@*/
void MPE_Thread_tls_destroy(MPE_Thread_tls_t * tls, int * err);

/*@
  MPE_Thread_tls_set - associate a value with the current thread in the thread local storage space
  
  Input Parameters:
+ tls - thread local storage space
- value - value to associate with current thread

  Output Parameter:
. err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_tls_set(MPE_Thread_tls_t * tls, void * value, int * err);

/*@
  MPE_Thread_tls_get - obtain the value associated with the current thread from the thread local storage space
  
  Input Parameter:
. tls - thread local storage space

  Output Parameters:
+ value - value associated with current thread
- err - location to store the error code; pointer may be NULL; error is zero for success, non-zero if a failure occurred
@*/
void MPE_Thread_tls_get(MPE_Thread_tls_t * tls, void ** value, int * err);


/*
 * Error values
 */
#define MPE_THREAD_SUCCESS MPE_THREAD_ERR_SUCCESS
#define MPE_THREAD_ERR_SUCCESS 0
/* FIXME: Define other error codes.  For now, any non-zero value is an error. */


/*
 * Implementation specific function definitions (usually in the form of macros)
 */
/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

/*
 * Threads
 */

/* MPE_Thread_create() defined in mpe_thread_posix.c */

#define MPE_Thread_exit()			\
{						\
    pthread_exit(NULL);				\
}

#define MPE_Thread_self(id_)			\
{						\
    *(id_) = pthread_self();			\
}

#define MPE_Thread_same(id1_, id2_, same_)			\
{								\
    *(same_) = pthread_equal(*(id1_), *(id2_)) ? TRUE : FALSE;	\
}

#define MPE_Thread_yield()						\
{									\
    /* FIXME: need to check for different types of yield */	\
    sched_yield();							\
}


/*
 *    Mutexes
 */

/* FIXME: using constant initializer if available */
#if !defined(MPICH_DEBUG_MUTEX) || !defined(PTHREAD_MUTEX_ERRORCHECK_NP)
#define MPE_Thread_mutex_create(mutex_ptr_, err_ptr_)                   \
{                                                                       \
    int err__;                                                          \
                                                                        \
    err__ = pthread_mutex_init((mutex_ptr_), NULL);                     \
    if ((err_ptr_) != NULL)                                             \
    {                                                                   \
	/* FIXME: convert error to an MPE_THREAD_ERR value */           \
	*(int *)(err_ptr_) = err__;                                     \
    }                                                                   \
}
#else /* MPICH_DEBUG_MUTEX */
#define MPE_Thread_mutex_create(mutex_ptr_, err_ptr_)                   \
{                                                                       \
    int err__;                                                          \
    pthread_mutexattr_t attr__;                                         \
                                                                        \
    pthread_mutexattr_init(&attr__);                                    \
    pthread_mutexattr_settype(&attr__, PTHREAD_MUTEX_ERRORCHECK_NP);    \
    err__ = pthread_mutex_init((mutex_ptr_), &attr__);                  \
    if (err__)                                                          \
        MPIU_Internal_sys_error_printf("pthread_mutex_init", err__,     \
                                       "    %s:%d\n", __FILE__, __LINE__);\
    if ((err_ptr_) != NULL)                                              \
    {                                                                    \
	/* FIXME: convert error to an MPE_THREAD_ERR value */            \
	*(int *)(err_ptr_) = err__;                                      \
    }                                                                    \
}
#endif

#define MPE_Thread_mutex_destroy(mutex_ptr_, err_ptr_)		\
{								\
    int err__;							\
								\
    err__ = pthread_mutex_destroy(mutex_ptr_);			\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#ifndef MPICH_DEBUG_MUTEX
#define MPE_Thread_mutex_lock(mutex_ptr_, err_ptr_)             \
{                                                               \
    int err__;                                                  \
    MPIU_DBG_MSG(THREAD,TYPICAL,"Enter MPE_Thread_mutex");      \
    err__ = pthread_mutex_lock(mutex_ptr_);                     \
    if ((err_ptr_) != NULL)                                     \
    {                                                           \
	/* FIXME: convert error to an MPE_THREAD_ERR value */   \
	*(int *)(err_ptr_) = err__;                             \
    }                                                           \
}
#else /* MPICH_DEBUG_MUTEX */
#define MPE_Thread_mutex_lock(mutex_ptr_, err_ptr_)             \
{                                                               \
    int err__;                                                  \
    MPIU_DBG_MSG(THREAD,TYPICAL,"Enter MPE_Thread_mutex");      \
    err__ = pthread_mutex_lock(mutex_ptr_);                     \
    if (err__)                                                  \
    {                                                           \
        MPIU_DBG_MSG_S(THREAD,TYPICAL,"  mutex lock error: %s", strerror(err__));       \
        MPIU_Internal_sys_error_printf("pthread_mutex_lock", err__,                     \
                                       "    %s:%d\n", __FILE__, __LINE__);              \
    }                                                          \
    if ((err_ptr_) != NULL)                                    \
    {                                                          \
	/* FIXME: convert error to an MPE_THREAD_ERR value */  \
	*(int *)(err_ptr_) = err__;                            \
    }                                                          \
}
#endif

#ifndef MPICH_DEBUG_MUTEX
#define MPE_Thread_mutex_unlock(mutex_ptr_, err_ptr_)           \
{                                                               \
    int err__;                                                  \
                                                                \
    MPIU_DBG_MSG(THREAD,TYPICAL,"Exiting MPE_Thread_mutex");    \
    err__ = pthread_mutex_unlock(mutex_ptr_);                   \
    if ((err_ptr_) != NULL)                                     \
    {                                                           \
	/* FIXME: convert error to an MPE_THREAD_ERR value */   \
	*(int *)(err_ptr_) = err__;                             \
    }                                                           \
}
#else /* MPICH_DEBUG_MUTEX */
#define MPE_Thread_mutex_unlock(mutex_ptr_, err_ptr_)           \
{                                                               \
    int err__;                                                  \
                                                                \
    MPIU_DBG_MSG(THREAD,TYPICAL,"Exiting MPE_Thread_mutex");    \
    err__ = pthread_mutex_unlock(mutex_ptr_);                   \
    if (err__)                                                  \
    {                                                           \
        MPIU_DBG_MSG_S(THREAD,TYPICAL,"  mutex unlock error: %s", strerror(err__));     \
        MPIU_Internal_sys_error_printf("pthread_mutex_unlock", err__,                   \
                                       "    %s:%d\n", __FILE__, __LINE__);              \
    }                                                           \
    if ((err_ptr_) != NULL)                                     \
    {                                                           \
	/* FIXME: convert error to an MPE_THREAD_ERR value */   \
	*(int *)(err_ptr_) = err__;                             \
    }                                                           \
}
#endif

#ifndef MPICH_DEBUG_MUTEX
#define MPE_Thread_mutex_trylock(mutex_ptr_, flag_ptr_, err_ptr_)    \
{                                                                    \
    int err__;                                                       \
                                                                     \
    err__ = pthread_mutex_trylock(mutex_ptr_);                       \
    *(flag_ptr_) = (err__ == 0) ? TRUE : FALSE;                      \
    if ((err_ptr_) != NULL)                                          \
    {                                                                \
	*(int *)(err_ptr_) = (err__ == EBUSY) : MPE_THREAD_SUCCESS ? err__;     \
	/* FIXME: convert error to an MPE_THREAD_ERR value */        \
    }                                                                \
}
#else /* MPICH_DEBUG_MUTEX */
#define MPE_Thread_mutex_trylock(mutex_ptr_, flag_ptr_, err_ptr_)    \
{                                                                    \
    int err__;                                                       \
                                                                     \
    err__ = pthread_mutex_trylock(mutex_ptr_);                       \
    if (err__ && err__ != EBUSY)                                     \
    {                                                                \
        MPIU_DBG_MSG_S(THREAD,TYPICAL,"  mutex trylock error: %s", strerror(err__));    \
        MPIU_Internal_sys_error_printf("pthread_mutex_trylock", err__,                  \
                                       "    %s:%d\n", __FILE__, __LINE__);              \
    }                                                                \
    *(flag_ptr_) = (err__ == 0) ? TRUE : FALSE;                      \
    if ((err_ptr_) != NULL)                                          \
    {                                                                \
	*(int *)(err_ptr_) = (err__ == EBUSY) : MPE_THREAD_SUCCESS ? err__; \
	/* FIXME: convert error to an MPE_THREAD_ERR value */        \
    }                                                                \
}
#endif

/*
 * Condition Variables
 */

#define MPE_Thread_cond_create(cond_ptr_, err_ptr_)		\
{								\
    int err__;							\
    								\
    err__ = pthread_cond_init((cond_ptr_), NULL);		\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#define MPE_Thread_cond_destroy(cond_ptr_, err_ptr_)		\
{								\
    int err__;							\
    								\
    err__ = pthread_cond_destroy(cond_ptr_);			\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#define MPE_Thread_cond_wait(cond_ptr_, mutex_ptr_, err_ptr_)		\
{									\
    int err__;								\
    									\
    /* The latest pthread specification says that cond_wait routines    \
       aren't allowed to return EINTR,	                                \
       but some of the older implementations still do. */		\
    do									\
    {									\
	err__ = pthread_cond_wait((cond_ptr_), (mutex_ptr_));		\
    }									\
    while (err__ == EINTR);						\
									\
    if ((err_ptr_) != NULL)						\
    {									\
	/* FIXME: convert error to an MPE_THREAD_ERR value */		\
	*(int *)(err_ptr_) = err__;					\
    }									\
}

#define MPE_Thread_cond_broadcast(cond_ptr_, err_ptr_)		\
{								\
    int err__;							\
    								\
    err__ = pthread_cond_broadcast(cond_ptr_);			\
								\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#define MPE_Thread_cond_signal(cond_ptr_, err_ptr_)		\
{								\
    int err__;							\
								\
    err__ = pthread_cond_signal(cond_ptr_);			\
								\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}


/*
 * Thread Local Storage
 */

#define MPE_Thread_tls_create(exit_func_ptr_, tls_ptr_, err_ptr_)	\
{									\
    int err__;								\
    									\
    err__ = pthread_key_create((tls_ptr_), (exit_func_ptr_));		\
									\
    if ((err_ptr_) != NULL)						\
    {									\
	/* FIXME: convert error to an MPE_THREAD_ERR value */		\
	*(int *)(err_ptr_) = err__;					\
    }									\
}

#define MPE_Thread_tls_destroy(tls_ptr_, err_ptr_)		\
{								\
    int err__;							\
    								\
    err__ = pthread_key_delete(*(tls_ptr_));			\
								\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#define MPE_Thread_tls_set(tls_ptr_, value_, err_ptr_)		\
{								\
    int err__;							\
								\
    err__ = pthread_setspecific(*(tls_ptr_), (value_));		\
								\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = err__;				\
    }								\
}

#define MPE_Thread_tls_get(tls_ptr_, value_ptr_, err_ptr_)	\
{								\
    *(value_ptr_) = pthread_getspecific(*(tls_ptr_));		\
								\
    if ((err_ptr_) != NULL)					\
    {								\
	/* FIXME: convert error to an MPE_THREAD_ERR value */	\
	*(int *)(err_ptr_) = MPE_THREAD_SUCCESS;		\
    }								\
}


#endif /* !defined(MPE_THREAD_H_INCLUDED) */
