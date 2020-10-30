def familiar_songs(unit_df):
    return unit_df.query("exposure == 'familiar' and call_type == 'song'")

def unfamiliar_songs(unit_df):
    return unit_df.query("exposure == 'unfamiliar' and call_type == 'song'")

def familiar_dcs(unit_df):
    return unit_df.query("exposure == 'familiar' and call_type == 'dc'")

def unfamiliar_dcs(unit_df):
    return unit_df.query("exposure == 'unfamiliar' and call_type == 'dc'")

def rewarded_songs(unit_df):
    return unit_df.query("reward_class == 'Rewarded' and call_type == 'song'")

def nonrewarded_songs(unit_df):
    return unit_df.query("reward_class == 'Nonrewarded' and call_type == 'song'")

def rewarded_dcs(unit_df):
    return unit_df.query("reward_class == 'Rewarded' and call_type == 'dc'")

def nonrewarded_dcs(unit_df):
    return unit_df.query("reward_class == 'Nonrewarded' and call_type == 'dc'")

def familiar_test_songs(unit_df):
    return unit_df.query("exposure == 'familiar' and relation == 'test' and call_type == 'song'")

def familiar_nontest_songs(unit_df):
    return unit_df.query("exposure == 'familiar' and relation != 'test' and call_type == 'song'")
    
def familiar_test_dcs(unit_df):
    return unit_df.query("exposure == 'familiar' and relation == 'test' and call_type == 'dc'")

def familiar_nontest_dcs(unit_df):
    return unit_df.query("exposure == 'familiar' and relation != 'test' and call_type == 'dc'")

def ripples(unit_df):
    return unit_df.query("call_type == 'ripple'")

def nonripples(unit_df):
    return unit_df.query("call_type != 'ripple'")

def dcs(unit_df):
    return unit_df.query("call_type == 'dc'")

def songs(unit_df):
    return unit_df.query("call_type == 'song'")
