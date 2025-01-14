using System;
using System.Text.Json;
using DisasterSimulation;

DirectoryInfo directoryInfo = new("../../../../source/faces");
IEnumerable<FileInfo> files = directoryInfo.EnumerateFiles("*", SearchOption.AllDirectories);
List<FaceData> faceData = new();
Console.WriteLine("面データ読み込み中...");
List<Task<FaceData[]>> readTasks = new();
foreach (FileInfo file in files)
{
    StreamReader streamReader = new(file.FullName);
    Task<FaceData[]> task = Task.Run(() =>
    {
        string json = streamReader.ReadToEnd();
        streamReader.Close();
        FaceData[] faceDataPart = JsonSerializer.Deserialize<FaceData[]>(json)!;
        return faceDataPart;
        //faceData.AddRange(faceDataPart);
    });
    readTasks.Add(task);
}
FaceData[][] faceDataParts = await Task.WhenAll(readTasks);
foreach (FaceData[] faceDataPart in faceDataParts)
{
    faceData.AddRange(faceDataPart);
}
faceData.Sort((a, b) => a.Index - b.Index);

DisasterSimulation.Simulator simulator = new(faceData.ToArray());
simulator.Start(70);
List<Result> result = simulator.result;

Console.WriteLine("シミュレーション結果を書き込み中...");
int lengthPerFile = 10;
List<Task> writeTasks = new();
List<String> fileNames = new();
for (int i = 0; i < result.Count; i += lengthPerFile)
{
    List<Result> resultPart = result.GetRange(i, Math.Min(lengthPerFile, result.Count - i));
    string resultJson = JsonSerializer.Serialize<Result[]>(resultPart.ToArray());
    StreamWriter streamWriter = new($"../../../../source/output/result{i}.json");
    fileNames.Add($"result{i}.json");
    Task task = Task.Run(() =>
    {
        streamWriter.Write(resultJson);
        streamWriter.Flush();
    });
    writeTasks.Add(task);
}
string fileNamesJson = JsonSerializer.Serialize(fileNames.ToArray());
StreamWriter fileNamesStreamWriter = new("../../../../source/output/fileNames.json");
Task fileNamesTask = Task.Run(() =>
{
    fileNamesStreamWriter.Write(fileNamesJson);
    fileNamesStreamWriter.Flush();
});
writeTasks.Add(fileNamesTask);
Task.WaitAll(writeTasks.ToArray());