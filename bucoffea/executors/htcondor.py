import htcondor

def htcondor_executor(items, function, accumulator, status=True, unit='items', desc='Processing',
                       function_args={}):
    for i, item in tqdm(enumerate(items), disable=not status, unit=unit, total=len(items), desc=desc):
        accumulator += function(item, **function_args)
    return accumulator
