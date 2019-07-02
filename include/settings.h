#ifndef _SETTINGS_H_
#define _SETTINGS_H_

char	sfn[80] = "./settings";

bool	read_settings (char *fn);
void	fprintf_settings (FILE *stream);

#endif	// _SETTINGS_H_

