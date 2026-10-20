param(
    [Parameter(Mandatory = $true)][string]$TestsDir,
    [Parameter(Mandatory = $true)][string]$ReportDir,
    [int]$BatchSize = 20,
    [int]$BatchTimeoutSeconds = 300
)

$ErrorActionPreference = 'Stop'
if ($BatchSize -le 0 -or $BatchTimeoutSeconds -le 0) {
    throw 'Batch size and timeout must be positive.'
}
$files = @(Get-ChildItem -LiteralPath $TestsDir -Recurse -File -Filter 'test_*.py' | Sort-Object FullName)
if ($files.Count -eq 0) {
    throw "No test_*.py files found in $TestsDir."
}
New-Item -ItemType Directory -Force -Path $ReportDir | Out-Null
$python = (Resolve-Path '.venv/Scripts/python.exe').Path
$batchCount = [int][Math]::Ceiling($files.Count / $BatchSize)

for ($index = 0; $index -lt $batchCount; $index++) {
    $first = $index * $BatchSize
    $last = [Math]::Min($first + $BatchSize - 1, $files.Count - 1)
    $batch = @($files[$first..$last])
    Write-Host "pytest batch $($index + 1)/$batchCount: $($batch.Count) files ($($batch[0].Name) through $($batch[-1].Name))"

    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = $python
    $startInfo.UseShellExecute = $false
    foreach ($argument in @('-m', 'pytest')) {
        [void]$startInfo.ArgumentList.Add($argument)
    }
    foreach ($file in $batch) {
        [void]$startInfo.ArgumentList.Add($file.FullName)
    }
    $report = Join-Path $ReportDir ('batch-{0:D2}.xml' -f ($index + 1))
    foreach ($argument in @('-v', '-x', '--timeout=90', '--timeout-method=thread', '--session-timeout=240', "--junitxml=$report")) {
        [void]$startInfo.ArgumentList.Add($argument)
    }
    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    if (-not $process.Start()) {
        throw "Could not start pytest batch $($index + 1)."
    }
    $start = [DateTime]::UtcNow
    while (-not $process.WaitForExit(30000)) {
        $elapsed = [int]([DateTime]::UtcNow - $start).TotalSeconds
        Write-Host "pytest batch $($index + 1) still running after ${elapsed}s (PID $($process.Id))"
        if ($elapsed -ge $BatchTimeoutSeconds) {
            try {
                $process.Kill($true)
                [void]$process.WaitForExit(10000)
            } catch {
                Write-Warning "Could not stop the pytest process tree: $_"
            }
            throw "pytest batch $($index + 1) exceeded ${BatchTimeoutSeconds}s."
        }
    }
    if ($process.ExitCode -ne 0) {
        throw "pytest batch $($index + 1) failed with exit code $($process.ExitCode)."
    }
}
Write-Host "All $($files.Count) Python binding test files passed."
