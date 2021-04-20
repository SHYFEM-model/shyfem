
void copy_record_to_array( int *ndim, int *ibuf, void *buf );

int get_next_record( int *ibuf, void **buf );
int get_next_record_32( int *ibuf, void **buf );
int get_next_record_64( int *ibuf, void **buf );
int skip_records( int nrec );
int copy_records( int nrec );

void set_read_file( char *string );
void set_write_file( char *string );
void close_read_file( void );
void close_write_file( void );
void init_read_file( void );

int *get_int( void *buf );
char *get_char( void *buf );
float *get_float( void *buf );

int read_record( void );
int read_record_32( void );
int read_record_64( void );
int read_record( void );
int write_record( void );
int write_record_32( void );
int write_record_64( void );

void set_swap_byte_order( int iflag );
void swap_byte_order( int n , void *buffer );
void swap_array( int n , int nbyte , void *buffer );

void skip_bytes( int nbytes );
void read_data( int nbytes );
void write_data( int nbytes );

void alloc_buffer( int nbytes );
void alloc_file_desc( void );

void Error( int ier, char *string );
void Error1( int ier, char *string, int i1 );
void Error2( int ier, char *string, int i1, int i2 );
void ErrorC( int ier, char *string, char *s );

