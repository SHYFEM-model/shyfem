#ifdef	DEBUG

#define	TRACE(format,args)		write(*,format) args 
#define	IFTRACE(condition,format,args)	if (condition) write(*,format) args 
#define	TESTNANVECTOR(args,s,g)		call testnan_vector(args,s,g)
#define	TESTNAN(scalar,s,g)		call testnan(scalar,s,g)
#define CHECKFLUX(A,B,C,D)              call unicflux(A,B,C,D)

#else

#define	TRACE(args)	
#define	IFTRACE(args)
#define	TESTNANVECTOR(args,s,g)
#define	TESTNAN(args,s,g)
#define CHECKFLUX(A,B,C,D)

#endif




