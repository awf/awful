function au_test_begin(tag)
% AU_TEST_BEGIN   Start a sequence of tests, and reset test stats

global all_test_status

fprintf(1, 'au_test_begin: start sequence [%s]\n', tag);

all_test_status.ok = 0;
all_test_status.failed = 0;

