#ifndef __GRAPHTYPES_H__
#define __GRAPHTYPES_H__

/* Read mapping struct */
#ifndef CONTIGNAME_LEN
#define CONTIGNAME_LEN  32
#endif

/* Chromosome position struct */
typedef struct {
    char        contigname[CONTIGNAME_LEN];
    uint64_t    contigstart;
    uint64_t    contigend;
    char        strand;
} chr_pos_t;

/* Position table entry */
typedef struct {
    chr_pos_t               pos;
    uint64_t                edgeindex;
} pos_entry_t;

/* Edge map entry */
typedef struct map_entry_s {
    pos_entry_t           * entry;
    struct map_entry_s    * next;
} map_entry_t;

#endif /* !defined(__DBTYPES_H__) */
