function au_test_result(ok)
% AU_TEST_RESULT   Report a test result

global all_test_status

if ok
    all_test_status.ok = all_test_status.ok + 1;
else
    all_test_status.failed = all_test_status.failed + 1;
end
