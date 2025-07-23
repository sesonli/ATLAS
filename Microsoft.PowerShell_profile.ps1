# Fix PSReadLine prediction issue that causes commands to hang
Set-PSReadLineOption -PredictionSource None

# Ensure we're using the latest PSReadLine version
Import-Module PSReadLine -RequiredVersion 2.3.6 -ErrorAction SilentlyContinue 