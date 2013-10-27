#ifndef _SAMUTIL_H
#define _SAMUTIL_H

int main_sam_sort(int argc, char ** argv);
int main_sam_checksort(int argc, char ** argv);
// int main_genome_to_transcript(int argc, char ** argv);
int main_tx2genome(int argc, char ** argv);
int main_get_tx_sequence(int argc, char ** argv);
int main_score_dist(int argc, char ** argv);
int main_score(int argc, char ** argv);
int main_sam_seqindex(int argc, char ** argv);
int main_sam_truncate(int argc, char ** argv);
int main_sam_rejoin(int argc, char ** argv);
int main_sam_filter(int argc, char ** argv);
int main_sam_extract(int argc, char ** argv);

//int main_generate_projection_header(int argc, char ** argv);

#endif // _SAMUTIL_H
