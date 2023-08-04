from qiime2.plugin import SemanticType

# Reg semantic types
PooledSampleData = SemanticType('PooledSampleData', field_names='type')

PooledSequencesWithQuality = SemanticType('PooledSequencesWithQuality',
                              variant_of=PooledSampleData.field['type'])

PooledSequences = SemanticType('PooledSequences',
                              variant_of=PooledSampleData.field['type'])

