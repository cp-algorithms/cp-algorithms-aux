const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const preview = require('./preview.cjs');

async function check({number = 7, state = 'open', fork = false, forkRun = false,
                      stale = false, link = false, git = false, symlinkRoot = false,
                      event = 'pull_request', conclusion = 'success'} = {}) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'preview-test-'));
  try {
    fs.writeFileSync(path.join(root, 'pr.json'), JSON.stringify({number}));
    fs.mkdirSync(path.join(root, 'site'));
    fs.writeFileSync(path.join(root, 'site/index.html'), '<p>preview</p>');
    if (link) fs.symlinkSync('../pr.json', path.join(root, 'site/link'));
    if (git) fs.mkdirSync(path.join(root, 'site/.git'));
    if (symlinkRoot) {
      fs.renameSync(path.join(root, 'site'), path.join(root, 'elsewhere'));
      fs.symlinkSync('elsewhere', path.join(root, 'site'));
    }
    const output = {};
    await preview({
      artifactRoot: root,
      core: {setOutput: (key, value) => { output[key] = value; }, info: () => {}},
      context: {repo: {owner: 'owner', repo: 'repo'}, payload: {workflow_run: {
        event, conclusion, head_sha: 'verified', head_repository: {full_name: forkRun ? 'fork/repo' : 'owner/repo'},
      }}},
      github: {rest: {pulls: {get: async ({pull_number}) => {
        assert.equal(pull_number, 7);
        return {data: {number: 7, state, base: {repo: {full_name: 'owner/repo'}},
          head: {sha: stale ? 'changed' : 'verified', repo: {full_name: fork ? 'fork/repo' : 'owner/repo'}}}};
      }}}},
    });
    return output;
  } finally {
    fs.rmSync(root, {recursive: true, force: true});
  }
}

(async () => {
  assert.deepEqual(await check(), {publish: 'true', number: '7', url: 'https://lib.cp-algorithms.com/7/'});
  for (const options of [{state: 'closed'}, {fork: true}, {forkRun: true}, {stale: true},
                         {event: 'push'}, {conclusion: 'failure'}]) {
    assert.deepEqual(await check(options), {publish: 'false'});
  }
  for (const number of ["7');throw new Error('injected');Number('7", '7$(id)', '../7', '7', -1, 1.5, 1e30]) {
    await assert.rejects(check({number}), /Invalid PR number/);
  }
  for (const options of [{link: true}, {git: true}, {symlinkRoot: true}]) {
    await assert.rejects(check(options), /Preview/);
  }
  console.log('17 preview identity and artifact checks passed');
})().catch(error => { console.error(error); process.exitCode = 1; });
