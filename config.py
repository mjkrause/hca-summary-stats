endpoints = {
    'hca_matrix_service_url': 'https://matrix.staging.data.humancellatlas.org/v1'
}

# Default parameters to requests an gene expression matrix.
post_expression_matrix_params = {
    'feature': 'gene',
    'fields': [
        'featurekey', 'featurename'
    ],
    'filter': {},
    'format': 'mtx'
}

s3_bucket = {'bucket_name': 'dev.project-assets.data.humancellatlas.org',
             'key': 'project-assets/project-stats/'}
