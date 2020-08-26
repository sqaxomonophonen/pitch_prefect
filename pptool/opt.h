#ifndef __OPT_H__

#ifdef OPT_STATIC
#define OPT_DEF static
#else
#define OPT_DEF extern
#endif

enum OptType {
	OPT_FLAG   = 1, // flags are booleans and take no arguments
	OPT_OPTION = 2, // options require an argument
	OPT_END = 0,
};
struct OptDef {
	enum OptType type;
	char short_opt;
	const char* long_opt;
	void* usr;
};

struct Opt {
	int is_invalid;
	int is_flag;
	int is_option;
	int is_switch;
	int is_npos;

	char* errmsg;

	char* arg;
	char* value;
	char short_opt;
	struct OptDef* def;
	int def_index;

	struct OptDef* defs;
	int argc;
	char** argv;
};

OPT_DEF void opt_init(struct Opt* opt, struct OptDef* defs, int argc, char** argv);
OPT_DEF int opt_next(struct Opt* opt);

#define __OPT_H__
#endif

#ifdef OPT_IMPLEMENTATION

#ifndef OPT_ASSERT
#include <assert.h>
#define OPT_ASSERT assert
#endif

#define OPT_UNREACHABLE OPT_ASSERT(!"UNREACHABLE")
#define OPT_INVALID_INPUT OPT_ASSERT(!"INVALID INPUT")

#include <string.h>

OPT_DEF void opt_init(struct Opt* opt, struct OptDef* defs, int argc, char** argv)
{
	memset(opt, 0, sizeof *opt);
	opt->defs = defs;
	opt->argc = argc;
	opt->argv = argv;
}

OPT_DEF int opt_next(struct Opt* opt)
{
	if (opt->argc == 0 || opt->is_invalid) return 0;

	int expect_value = 0;

	opt->def_index = -1;
	opt->def = NULL;
	opt->short_opt = 0;
again:

	opt->is_invalid = 0;
	opt->is_flag = 0;
	opt->is_option = 0;
	opt->is_switch = 0;
	opt->is_npos = 0;
	opt->errmsg = NULL;
	opt->value = NULL;

	if (expect_value && opt->argc == 0) {
		opt->is_invalid = 1;
		opt->errmsg = "expected value for option";
		return 1;
	} else {
		OPT_ASSERT(opt->argc > 0);
	}

	char* arg = *opt->argv;
	opt->arg = arg;
	opt->argc--;
	opt->argv++;
	int len = strlen(arg);
	int n_dashes;
	for (n_dashes = 0; arg[n_dashes] == '-'; n_dashes++);

	int do_search_def = 0;
	char short_opt = 0;
	char* long_opt = NULL;
	size_t long_opt_sz = 0;

	char* value = NULL;

	if (n_dashes == 1 && len > 1) {
		if (len == 2) {
			short_opt = arg[1];
			do_search_def = 1;
		} else {
			opt->is_invalid = 1;
			opt->errmsg = "invalid short opt";
			return 1;
		}
	} else if (n_dashes == 2 && len > 2) {
		int argend;
		for (argend = 2; argend < len && arg[argend] != '='; argend++);
		long_opt = &arg[2];
		long_opt_sz = argend - 2;
		if (argend != len) value = &arg[argend+1];
		do_search_def = 1;
	} else {
		value = arg;
	}

	if (expect_value) {
		if (!value) {
			opt->is_invalid = 1;
			opt->errmsg = "expected argument value";
			return 1;
		}
		opt->is_option = 1;
		opt->is_switch = 1;
		opt->value = value;
		return 1;
	} else if (do_search_def) {
		int found_def_at_index = -1;
		struct OptDef* def;
		int i = 0;
		for (def = opt->defs; def->type; def++, i++) {
			if (short_opt) {
				if (short_opt == def->short_opt) {
					found_def_at_index = i;
					break;
				}
			} else if (long_opt) {
				if (strlen(def->long_opt) == long_opt_sz && memcmp(long_opt, def->long_opt, long_opt_sz) == 0) {
					found_def_at_index = i;
					break;
				}
			} else {
				OPT_UNREACHABLE;
			}
		}

		if (found_def_at_index == -1) {
			opt->is_invalid = 1;
			opt->errmsg = "flag/option not defined";
			return 1;
		}

		OPT_ASSERT(def != NULL);

		opt->def_index = found_def_at_index;
		opt->def = def;
		opt->short_opt = def->short_opt;

		if (def->type == OPT_FLAG) {
			opt->is_flag = 1;
			opt->is_switch = 1;
			return 1;
		} else if (def->type == OPT_OPTION) {
			if (value == NULL) {
				expect_value = 1;
				goto again;
			} else {
				opt->is_option = 1;
				opt->is_switch = 1;
				opt->value = value;
				return 1;
			}
		} else {
			OPT_INVALID_INPUT;
		}

		return 1;
	} else if (value) {
		opt->is_npos = 1;
		opt->value = value;
		return 1;
	}

	OPT_UNREACHABLE;
}

#endif

#ifdef OPT_TEST

// run unit test:
//  $ cc -x c -std=c99 -O2 -Wall -DOPT_IMPLEMENTATION -DOPT_TEST opt.h -o test_opt && ./test_opt
// or use GDB:
//  $ cc -x c -std=c99 -O0 -g -Wall -DOPT_IMPLEMENTATION -DOPT_TEST opt.h -o test_opt && gdb ./test_opt

#define ARRAY_LENGTH(x) (sizeof(x) / sizeof(x[0]))

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int _argc, char** _argv)
{
	{
		// FIXME not much of a test currently, more like a sample program :-)
		char* argv[] = {"skip me", "npos1", "npos2", "-y", "npos3", "-x", "npos4", "-z", "--longz", "--longy", "--long-x", "--sample-rate", "48000"};
		//char* argv[] = {"skip me", "npos1", "npos2", "-y", "npos3", "-x", "npos4", "-z", "--longz", "--longy", "--long-x", "--sample-rate"};
		int argc = ARRAY_LENGTH(argv);

		struct OptDef defs[] = {
			{OPT_FLAG, 'x', "long-x"},
			{OPT_FLAG, 'y', "longy"},
			{OPT_FLAG, 'z', "longz"},
			{OPT_OPTION, 'S', "sample-rate"},
			{0},
		};
		struct Opt opt;
		opt_init(&opt, defs, argc-1, argv+1);
		while (opt_next(&opt)) {
			assert(!opt.is_invalid);
			if (opt.is_flag) {
				printf("flag index=%d short=%c\n", opt.def_index, opt.short_opt);
			} else if (opt.is_option) {
				printf("option index=%d short=%c value=%s\n", opt.def_index, opt.short_opt, opt.value);
			} else if (opt.is_npos) {
				printf("npos value=%s\n", opt.value);
			} else {
				OPT_UNREACHABLE;
			}
		}
	}

	return EXIT_SUCCESS;
}

#endif
