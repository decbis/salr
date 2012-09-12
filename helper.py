from salr import character_count_tester


def _read(token_type):
    lst = []

    f = open('./tokens.txt', 'r')
    for line in f:
        line_lst = line.split('"')
        if line_lst[1] == token_type:
            lst.append(line_lst[3])
    f.close()

    return lst


def _chart(histogram):
    for k in sorted(histogram.keys()):
        v = histogram[k]
        print k, ' ', '*' * v


nonces = _read('session_nonce')
auth_tokens = _read('auth_token')

test_results = character_count_tester(nonces)
#character_count_tester(auth_tokens)

for (prob, histogram) in test_results:
    print prob
    _chart(histogram)
