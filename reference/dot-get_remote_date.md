# Get Remote File Modification Date

Retrieves the last modification date of a remote file from its HTTP
headers. Useful for checking if cached data needs updating.

## Usage

``` r
.get_remote_date(url)
```

## Arguments

- url:

  Character string specifying the URL to check. If NULL, returns NULL.

## Value

A Date object representing the last modification time, or NA if the
request fails. Returns NA if the URL is NULL or if the remote server is
unreachable.

## Examples

``` r
if (FALSE) { # \dontrun{
  .get_remote_date("https://example.com/data.txt")
} # }
```
