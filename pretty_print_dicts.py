

constants_dict = {
    'Oregon': 'OR',
    'Florida': 'FL',
    'California': 'CA',
    'New York': 'NY',
    'Michigan': 'MI'
}

uniform_len = 50
spaces = " "*uniform_len

for param in constants_dict:
    padded_param = param[:uniform_len]+spaces[:uniform_len-len(param)]
    param_val = str(constants_dict[param])
    print padded_param+":"+param_val
